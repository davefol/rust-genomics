#![crate_name = "genomics"]
use ndarray;
use ndarray::ShapeBuilder;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Mutex};

pub mod prelude;

pub mod observable;
pub mod index_of_association;

pub type Groups = HashMap<String, Arc<Group>>;
pub type Meta = HashMap<String, String>;
pub type AlleleCount = u32;
pub type Variations = Arc<Mutex<BTreeMap<String, Arc<Variation>>>>;
pub type Loci = BTreeMap<String, Arc<Locus>>;
pub type Allele = (Arc<Locus>, Arc<Variation>);
pub type Individuals = BTreeMap<String, Individual>;
pub type Genome = HashMap<Allele, AlleleCount>;

pub trait LociExt {
    fn n_alleles(&self) -> usize;
}

impl LociExt for Loci {
    fn n_alleles(&self) -> usize {
        self.iter().fold(0, |count, (_, locus)| {
            count + locus.variations.lock().unwrap().len()
        })
    }
}

#[derive(Hash, Eq, PartialEq)]
pub struct Variation {
    name: String,
}

impl Variation {
    pub fn new(name: &str) -> Self {
        Self { name: name.into() }
    }
}

pub enum LocusHint {
    Classical,
    Microsatellite,
}

pub struct Locus {
    name: String,
    variations: Variations,
    hint: LocusHint,
}

impl Hash for Locus {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}

impl PartialEq for Locus {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl Eq for Locus {}

impl Locus {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.into(),
            variations: Variations::new(Mutex::new(BTreeMap::new())),
            hint: LocusHint::Microsatellite,
        }
    }
}

#[derive(Hash, PartialEq, Eq)]
pub struct Group {
    name: String,
}

impl Group {
    pub fn new(name: &str) -> Self {
        Self { name: name.into() }
    }
}

pub struct Individual {
    name: String,
    genome: Genome,
    groups: HashSet<Arc<Group>>,
    meta: Meta,
}

impl Individual {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.into(),
            genome: Genome::new(),
            groups: HashSet::new(),
            meta: Meta::new(),
        }
    }
}

pub struct AlleleMatrix {
    data: ndarray::Array2<AlleleCount>,
    loci: Vec<(usize, usize)>,
    dirty: bool,
}

impl AlleleMatrix {
    pub fn new() -> Self {
        Self {
            data: ndarray::Array2::<AlleleCount>::zeros((0, 0)),
            loci: vec![],
            dirty: false,
        }
    }

    pub fn from_vec(individuals: usize, loci: Vec<(usize, usize)>, data: Vec<AlleleCount>) -> Result<Self, Box<dyn Error>> {
        let alleles = data.len() / individuals;
        Ok(Self {
            data: ndarray::Array::from_shape_vec((individuals, alleles).strides((alleles, 1)), data)?,
            loci: loci,
            dirty: false,
        })
    }

    /// Computes the frequency matrix from allele counts
    pub fn frequency(&self) -> Result<ndarray::Array2<f32>, Box<dyn Error>> {
        let mut freqs = ndarray::Array2::from_elem(self.data.dim(), 0.0);
        let loci: [(usize, usize); 2] = [(0,3), (3, 6)];
    
        ndarray::Zip::from(freqs.genrows_mut())
        .and(self.data.genrows())
        .apply(|mut freqs, row| {
            let x = ndarray::Array::from(loci.iter().flat_map(|(start, end)| {
                let loci_sum = row.slice(ndarray::s![*start..*end]).sum() as f32;
                row.slice(ndarray::s![*start..*end]).map(|x| *x as f32 / loci_sum).to_vec()
            }).collect::<Vec<_>>());
            freqs.assign(&x);
        });
        Ok(freqs)
    }
}

/// An observation of an Individual
pub enum Observation {
    /// An Observation that an individual has an Allele
    /// Individual's name, Locus's name, Variation's name
    Allele(String, String, String),

    /// An Observation that an Individual belongs to a Group
    /// Individual's name, Group's name
    Group(String, String),

    /// An `Observation` that an `Individual` has associated metadata
    /// Individual's name, Meta data description, Meta data content. 
    Meta(String, String, String),
}

pub struct Sample {
    loci: Loci,
    groups: Groups,
    individuals: Individuals,
    matrix: AlleleMatrix,
}

impl Sample {
    /// Constructs a new empty `Sample`
    ///
    /// The `Sample` can be filled up iteratively by calling
    /// `observation()`.
    pub fn new() -> Self {
        Self {
            loci: Loci::new(),
            groups: Groups::new(),
            individuals: Individuals::new(),
            matrix: AlleleMatrix::new(),
        }
    }

    /// Recreates the `Sample`'s `matrix`.
    ///
    /// This function is called before a matrix calculation
    /// so there is no need to explicitly call it after observing data.
    pub fn flush(&mut self) -> Result<(), Box<dyn Error>> {

        let loci: Vec<(usize, usize)> = self.loci.iter().enumerate().map(|(i, x)| {
            (i, i + x.1.variations.lock().unwrap().len())
        }).collect();
        self.matrix = AlleleMatrix::from_vec(self.individuals.len(), loci, Vec::<AlleleCount>::from(&*self))?;
        Ok(())
    }

    /// Returns the allele in this locus
    ///
    /// This will create the locus and allele if needed and mark
    /// the Sample as dirty. Dirty samples will updated before
    /// the next call that requires a matrix calculation or when
    /// `flush()` is called.
    pub fn allele(&mut self, locus: &str, variation: &str) -> Allele {
        (
            self.loci
                .entry(locus.into())
                .or_insert({
                    self.matrix.dirty = true;
                    Arc::new(Locus::new(locus))
                })
                .clone(),
            self.loci
                .get(locus)
                .unwrap() // Guaranteed to exist since we just inserted it.
                .variations
                .lock()
                .unwrap()
                .entry(variation.into())
                .or_insert({
                    self.matrix.dirty = true;
                    Arc::new(Variation::new(variation))
                })
                .clone(),
        )
    }

    /// Returns a reference to a `Group`
    ///
    /// This will create the `Group` if needed and mark the `Sample`
    /// as dirty.
    pub fn group(&mut self, group: &str) -> Arc<Group> {
        self.groups
            .entry(group.into())
            .or_insert({
                self.matrix.dirty = true;
                Arc::new(Group::new(group.into()))
            })
            .clone()
    }

    /// Observes a single `Observation`
    ///
    /// This function is normally called by a `Sample`'s observe
    /// function to read in data. To read in data from an arbitrary
    /// data source. Implement an Iterator with type Item = Observation.
    pub fn _observe(&mut self, observation: Observation) {
        match &observation {
            Observation::Allele(individual, locus, variation) => {
                let allele = self.allele(&locus, &variation);
                *self
                    .individuals
                    .entry(individual.into())
                    .or_insert(Individual::new(individual))
                    .genome
                    .entry(allele)
                    .or_insert(0) += 1;
            }
            Observation::Group(individual, group) => {
                let arc_group = self.group(group);
                self.individuals
                    .entry(individual.into())
                    .or_insert(Individual::new(individual))
                    .groups
                    .insert(arc_group);
            },
            Observation::Meta(individual, meta, content) => {
                self.individuals
                    .entry(individual.into())
                    .or_insert(Individual::new(individual))
                    .meta
                    .insert(meta.into(), content.into());
            }
        }
    }

    /// Observe all the data in the argument.
    pub fn observe<I>(&mut self, observable: I) -> Result<(), Box<dyn Error>>
    where 
        I: Iterator<Item = Result<Observation, Box<dyn Error>>> 
    {
        for observation in observable {
            self._observe(observation?);
        }
        Ok(())
    }

    /// A list of the names of all loci in a sample
    pub fn loci_names(&self) -> Vec<&String> {
        self.loci.keys().collect()
    }

    pub fn variations(&self, locus: &str) -> Option<Vec<String>> {
        self.loci.get(locus).map(|loc| {
            loc.variations.lock().unwrap().keys().map(|x| x.to_string()).collect()
        })
    }
}

impl From<&Sample> for Vec<AlleleCount> {
    fn from(sample: &Sample) -> Vec<AlleleCount> {
        let mut v = vec![];
        for (_, individual) in sample.individuals.iter() {
            for (_, locus) in sample.loci.iter() {
                for (_, variation) in locus.variations.lock().unwrap().iter() {
                    v.push(
                        match individual.genome.get(&(locus.clone(), variation.clone())) {
                            None => 0,
                            Some(c) => *c,
                        },
                    )
                }
            }
        }
        v
    }
}

