Protein Intensity Mapper (PrIntMap-R)
=====================

# Overview
PrIntMap-R is a suite of visualization tools designed to help with the analysis of shotgun proteomic data at the protein coverage level. The plots generated here can help to visualize uneven coverage of a protein, compare two or multiple samples, and look at the different peptides that contribute to the overall coverage of the protein. For use examples, as well as example data for download, check out the `Examples` tab. To use the tools, click on the `Run` tab.  

-------

# Inputs

### Database File
The database file should contain all the proteins that you many be interested in, and generally should be the same file used in your database search. This file should be in [UniProt fasta format](https://www.uniprot.org/help/fasta-headers) which looks like this:
```
>db|UniqueIdentifier|EntryName ProteinName OS=OrgName OX=OrgID [GN=GGN ]PE=PE SV=SV
AMINOACIDSEQUENCEPEPTIDE...
```
The program will ask for an `Accession ID`. This generally is the `UniqueIdentifier` in the example above, but it will actually match anything in the `UniqueIdentifier` or `EntryName` fields.

### Peptide Files
The peptide files are exported from the database search. These are either `.csv`, `.tsv`, or `.txt` files. Currently supported formats are from [PEAKS](https://www.bioinfor.com/peaks-online/), [MSFragger](https://msfragger.nesvilab.org/), and [MaxQuant](https://www.maxquant.org/), with plans to include Proteome Discover and MetaMorpheus soon.  

* PEAKS:  
* MSFragger:  
* MaxQuant:
* Proteome Discover:
  * *Coming Soon*
* MetaMorpheus:
  * *Coming Soon*

### Other Options

-------

# Outputs

### Basic

### Two Samples

### Multiple Samples

### Annotation

### Unique Peptides

### Stacked Peptide Plot

### Percent Coverage

### Export

------

<img src="Champion_lab_v1.png" alt="drawing" width="130"/>  
[Champion Lab]("https://championlab.weebly.com/)  
University of Notre Dame  

