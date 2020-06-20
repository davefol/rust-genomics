use ndarray;

pub struct Variation {
    id: u16,
    name: String,
}

pub struct Locus {
    id: u16,
    name: String,
    variations: Vec<Variation>,
}

pub struct Allele<'locus> {
    locus: &'locus Locus,
    variation: &'locus Variation,
    count: u8,
}

pub struct Group {
    id: u32,
    name: String,
}

pub struct Individual<'locus, 'group> {
    id: u32,
    name: Option<String>,
    alleles: Vec<Allele<'locus>>,
    group: Option<&'group Group>,
}

pub struct Sample<'locus, 'group> {
    individuals: Vec<Individual<'locus, 'group>>,
    array: ndarray::Array2<u8>,
}
