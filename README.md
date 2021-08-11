## Molecules standardization using CGRtools package

This app uses [CGRtools package](https://github.com/stsouko/CGRtools) to standardize molecular data. 
It takes *SDF* or *SMILES* input files, checks their structure on valence and aromatic ring errors and
applies standardization rules in the folowing order:

1. Convertation of molecule to the kekule form
2. Standardization of functional groups
3. Hidding of explicit hydrogens
4. Convertation of molecule to the aromatic (thiele) form

For additional information, 
please, refer to the [CGRtools documentation](https://cgrtools.readthedocs.io/tutorial/3_standardization.html).