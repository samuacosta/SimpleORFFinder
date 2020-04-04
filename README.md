# SimpleORFFinder

A simple DNA open reading frame finder with custom genetic code support.

## Usage

```console
python simple_orf_finder.py [-h] [-i INPUT [INPUT ...]] [-o OUTPUT] [-min MINIMUM] [-g GENCODE] [-start CODON [CODON ...]] [-stop CODON [CODON ...]] [-f FRAME [FRAME ...]]
```

Finds open reading frames from DNA nucleotide sequences provided using INPUT file(s) in FASTA format, with MINIMUM nucleotide length, based on a GENCODE CSV file, using provided start/stop CODON(S) in FRAME(S) reading frames. Output is also in FASTA format.

## Features

The script searches for all the possible open reading frames from DNA nucleotide sequences and writes them as a FASTA formatted file, to the standard output or to a file if specified. It supports multiple files and multiple sequences per file. Custom genetic code can be supplied via a CSV file (codon, amino acid) and also start/stop codon(s). Minimum ORF length and limitation to specific subset of reading frames also supported.

In case of nested ORFs, the longest is used. Input must be in FASTA format. Unusual codons are supported, the resulting unknown amino acid is identified as X and a warning is shown in standard output. ORFs containing 	ambiguous nucleotides are ignored.

ORFs are labelled combining the start of the name provided in the input FASTA header with the number of frame and number of each ORF within the frame. Length, 	position and frame are also written in the FASTA header.

## Examples
Minimal usage (test_input.fasta as input, 150 minimum nucleotide length, standard genetic code and start/stop codons, all frames, prints to stdout):
```console
python simple_orf_finder.py
```

Full parameter specification (input/output file(s), minimum ORF length of 100, custom genetic code file, custom start/stop codons, specifies the three forward frames and the three backwards frames to be considered in the search):
```console
python simple_orf_finder.py -i my_first_input.fasta my_second_input.fasta -o my_output.fasta -min 100 -g custom_genetic_code.csv -start ATG -stop TAG TAA TGA -f 1 2 3 4 5 6
```

Print command usage:
```console
python simple_orf_finder.py -h
```

## Options
```console
-i [INPUT [INPUT ...]], --input INPUT [INPUT ...]
	The name of the input file(s) in FASTA format.
	[Default=test_input.fasta]
-o OUTPUT, --output OUTPUT
	The name of the optional output file.
	[Default=stdout]
-min MINIMUM, --minimum MINIMUM
	The minimum size of ORFs to search, in nucleotides.
	[Default=150]
-g GENCODE, --gencode GENCODE
	The name of the file with the genetic code to use.
	[Default=standard_genetic_code.csv]
-start CODON, --start CODON
	List of possible start codons to use.
	[Default=ATG]
-stop CODON, --stop CODON
	List of possible stop codons to use.
	[Default=TAG,TAA,TGA]
-f [FRAMES [FRAMES ...]], --frames [FRAMES [FRAMES ...]]
	The frames that should be considered (1 to 6).
	[Default=All]
-h, --help
	Shows a brief “usage” text of the command and exits
```

## Sample output
```console
>I.CLAUDIUS_F1_0001 FRAME1 203 1888
MQVSRRKFFKICAGGMAGTSAAMLGFAPANVLAAPREYKLLRAFESRNTCTYCAVSCGMLLYSTGKPYNS
LSSHTGTNTRSKLFHIEGDPDHPVSRGALCPKGAGSLDYVNSESRSLYPQYRAPGSDKWERISWKDAIKR
IARLMKDDRDANFVEKDSNGKTVNRWATTGIMTASAMSNEAALLTQKWIRMLGMVPVCNQANT
```

## Author
Samuel Acosta [samuel.acostamelgarejo AT postgrad.manchester.ac.uk]

## License
[MIT](https://choosealicense.com/licenses/mit/)