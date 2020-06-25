use crate::prelude::*;
use csv;
use std::collections::{HashSet, VecDeque};
use std::error::Error;
use std::io::Read;

enum ObservationPartial {
    Allele(String, String),
    Group(String),
    Meta(String, String),
}

impl ObservationPartial {
    fn to_observation(&self, individual: &str) -> Observation {
        match self {
            Self::Allele(locus, variation) => Observation::Allele(
                individual.to_string(),
                locus.to_string(),
                variation.to_string(),
            ),
            Self::Group(group) => Observation::Group(individual.to_string(), group.to_string()),
            Self::Meta(meta, content) => Observation::Meta(
                individual.to_string(),
                meta.to_string(),
                content.to_string(),
            ),
        }
    }
}

#[derive(Clone)]
enum Field {
    Locus(String),
    Name,
    Group,
    GroupPresence(String),
    Meta(String),
}

/// Produces Observations from u8 delimted data
///
/// `Csv` implements Iterator so it can be passed
/// directly to `Sample::observe()`
pub struct Csv {
    records: std::iter::Enumerate<csv::StringRecordsIntoIter<Box<dyn Read>>>,
    fields: Option<Vec<Field>>,
    separator: String,
    observation_buffer: VecDeque<Observation>,
    group_presence_identifier: String,
}

impl Csv {
    fn new(
        records: csv::StringRecordsIntoIter<Box<dyn Read>>,
        fields: Option<Vec<Field>>,
        separator: &str,
        group_presence_identifier: &str,
    ) -> Self {
        Self {
            records: records.into_iter().enumerate(),
            fields,
            separator: separator.to_owned(),
            observation_buffer: VecDeque::new(),
            group_presence_identifier: group_presence_identifier.to_owned(),
        }
    }
}

impl Iterator for Csv {
    type Item = Result<Observation, Box<dyn Error>>;

    fn next(&mut self) -> Option<Result<Observation, Box<dyn Error>>> {
        if self.observation_buffer.len() == 0 {
            match self.records.next() {
                None => {
                    return None;
                }
                Some((idx, Ok(row))) => {
                    let mut individual = idx.to_string();
                    let mut partials = vec![];
                    if let Some(fields) = &self.fields {
                        for (i, field) in row.iter().enumerate() {
                            match &fields[i] {
                                Field::Name => {
                                    individual = field.to_string();
                                }
                                Field::Locus(s) => {
                                    for x in field.split(&self.separator) {
                                        partials
                                            .push(ObservationPartial::Allele(s.into(), x.into()));
                                    }
                                }
                                Field::Group => {
                                    partials.push(ObservationPartial::Group(field.into()));
                                }
                                Field::GroupPresence(s) => {
                                    if field == self.group_presence_identifier {
                                        partials.push(ObservationPartial::Group(s.into()));
                                    }
                                }
                                Field::Meta(s) => {
                                    partials.push(ObservationPartial::Meta(s.into(), field.into()));
                                }
                            }
                        }
                    } else {
                        for (i, field) in row.iter().enumerate() {
                            for x in field.split(&self.separator) {
                                partials.push(ObservationPartial::Allele(i.to_string(), x.into()));
                            }
                        }
                    }
                    self.observation_buffer = partials
                        .iter()
                        .map(|x| x.to_observation(&individual))
                        .collect();
                }
                Some((idx, Err(_))) => {
                    return None;
                }
            }
        }

        Some(Ok(self.observation_buffer.pop_front().unwrap()))
    }
}

pub struct CsvBuilder {
    headers: bool,
    delimiter: u8,
    separator: String,
    name_field: Option<String>,
    group_fields: HashSet<String>,
    group_field: Option<String>,
    meta_fields: HashSet<String>,
    group_presence_identifier: String,
}

impl CsvBuilder {
    /// Construct a new Csv builder
    pub fn new() -> Self {
        Self {
            headers: true,
            delimiter: b',',
            separator: "/".to_owned(),
            name_field: None,
            group_fields: HashSet::new(),
            group_field: None,
            meta_fields: HashSet::new(),
            group_presence_identifier: "Y".to_owned(),
        }
    }

    pub fn headers(&mut self, headers: bool) -> &mut Self {
        self.headers = headers;
        self
    }

    pub fn delimiter(&mut self, delimiter: u8) -> &mut Self {
        self.delimiter = delimiter;
        self
    }

    pub fn separator(&mut self, separator: &str) -> &mut Self {
        self.separator = separator.to_owned();
        self
    }

    pub fn name_field(&mut self, name_field: &str) -> &mut Self {
        self.name_field = Some(name_field.to_owned());
        self
    }

    pub fn group_field(&mut self, group_field: &str) -> &mut Self {
        self.group_field = Some(group_field.to_owned());
        self
    }

    pub fn group_fields(
        &mut self,
        group_fields: HashSet<String>,
        group_presence_identifier: &str,
    ) -> &mut Self {
        self.group_fields = group_fields;
        self.group_presence_identifier = group_presence_identifier.to_owned();
        self
    }

    pub fn meta_fields(&mut self, meta_fields: HashSet<String>) -> &mut Self {
        self.meta_fields = meta_fields;
        self
    }

    pub fn from_reader(&self, reader: Box<dyn Read>) -> Result<Csv, Box<dyn Error>> {
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(self.headers)
            .delimiter(self.delimiter)
            .from_reader(reader);

        let fields = if self.headers {
            Some(
                rdr.headers()?
                    .iter()
                    .map(|s| {
                        if let Some(name_field) = &self.name_field {
                            if s == name_field {
                                return Field::Name;
                            }
                        }

                        if self.group_fields.contains(s) {
                            return Field::GroupPresence(s.into());
                        }

                        if let Some(group_field) = &self.group_field {
                            if group_field == s {
                                return Field::Group;
                            }
                        }

                        if self.meta_fields.contains(s) {
                            return Field::Meta(s.into());
                        }

                        Field::Locus(s.into())
                    })
                    .collect(),
            )
        } else {
            None
        };

        Ok(Csv::new(
            rdr.into_records(),
            fields,
            &self.separator,
            &self.group_presence_identifier,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn test_csv_with_header_has_correct_loci() -> Result<(), Box<dyn Error>> {
        let mut sample = Sample::new();
        sample.observe(CsvBuilder::new().from_reader(Box::new("test\n0/0".as_bytes()))?);
        assert_eq!(sample.loci_names(), vec!["test"]);
        Ok(())
    }

    #[test]
    fn test_csv_has_correct_variations() -> Result<(), Box<dyn Error>> {
        let mut sample = Sample::new();
        sample.observe(
            CsvBuilder::new().from_reader(Box::new("test\n0/1/2/3/4\n5/6/7/8/9".as_bytes()))?,
        );
        assert_eq!(
            sample.variations("test").unwrap(),
            (0..=9).map(|x| x.to_string()).collect::<Vec<String>>()
        );
        Ok(())
    }
}
