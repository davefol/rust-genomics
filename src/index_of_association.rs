use crate::prelude::*;
use ndarray;
use std::collections::HashMap;
use std::error::Error;

pub struct IndexOfAssociationSummary {
    index_of_association: f32,
}

pub trait IndexOfAssociation {
    fn index_of_association(&mut self) -> Result<IndexOfAssociationSummary, Box<dyn Error>>;
}

impl IndexOfAssociation for Sample {
    fn index_of_association(&mut self) -> Result<IndexOfAssociationSummary, Box<dyn Error>> {
        if self.matrix.dirty {
            self.flush()?;
        }

        let freqs = self.matrix.frequency()?;
        let n_freqs = freqs.shape()[0];

        let mut indices: HashMap<(usize, usize), usize> = HashMap::new();
        let mut counter = 0;
        for i in 0..n_freqs - 1 {
            for j in i..n_freqs {
                indices.insert((i, j), counter);
                counter += 1;
            }
        }

        let n_distances = n_freqs * (n_freqs - 1) / 2;
        let n_loci = self.matrix.loci.len();
        let mut distances = ndarray::Array::zeros((n_distances, n_loci));

        (0..n_freqs - 1).for_each(|i| {
            ((i + 1)..n_freqs).for_each(|j| {
                self.matrix
                    .loci
                    .iter()
                    .enumerate()
                    .for_each(|(idx, (start, end))| {
                        distances[[indices[&(i, j)], idx]] =
                            (&freqs.row(i).slice(ndarray::s![*start..*end])
                                - &freqs.row(j).slice(ndarray::s![*start..*end]))
                                .map(|x| x.abs())
                                .sum();
                    });
            })
        });

        let variance = (ndarray::Zip::from(distances.genrows())
            .apply_collect(|row| row.sum().powf(2.0))
            .sum()
            - ndarray::Zip::from(distances.genrows())
                .apply_collect(|row| row.sum())
                .sum()
                .powf(2.0)
                / n_distances as f32)
            / n_distances as f32;

        let expected_variance: f32 = (0..n_loci)
            .map(|n| {
                (distances.column(n).map(|x| x.powf(2.0)).sum()
                    - (distances.column(n).sum() / n_distances as f32))
                    / n_distances as f32
            })
            .sum();

        let index_of_association =  (variance / expected_variance) - 1.0;

        Ok(IndexOfAssociationSummary {
            index_of_association: index_of_association,
        })
    }
}

#[cfg(test)]
mod tests {

    //#[test]
    //fn test_index_of_association() -> Result<(), Box<dyn Error>> {
        
    //}
}
