# Construction Phylogenetic Trees of ALDH2 Among Species

## Dependencies

+ `clustalw2`, the latest version of which can be [downloaded here](http://www.clustal.org/download/current/)
    + Make sure this is on your system path for ease of use, especially on windows and linux.
+ Python dependencies, which can be installed with `pip install -r requirements.txt`. See the [pip documentation](https://pip.pypa.io/en/stable/quickstart/) for more assistance.

## Usage

+ Retrieve Amino Acid sequences to use.
+ Place the AA sequences, in FASTA format, in `sequences_aa.fa`
+ Run `python3 team2-final-project.py` with the following optional flags at the end.
    + `-h` prints usage help text and quits
    + `-n` skips `clustalw2` alignment and instead uses an existing `sequences_aa.aln` file