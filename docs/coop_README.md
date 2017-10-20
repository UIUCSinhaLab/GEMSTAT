# Cooperativities

The user can provide a file specifying TF-TF cooperative binding interactions (including antagonism) with the

-c coop_file.txt

command line argument. The file can have any name and extension.

The basic format of the file is as follows

```
onetf	anothertf
yetanother	andanother
selfinteracing	selfinteracting
```

Where each line has a tab separated name of two TFs that interact. (Self interaction is handled by naming the same TF twice.)


## Special cooperativities or distance limits.

All cooperativities will be subject to the global cutoff specified by the maximum of the -ct and the -rt options. (For algorithmic reasons.)

Special cooperativities (differning from the global default) are specified in the coop.txt file.

```
first   first   DIMER   11      +       -
```

Asks for the tf named `first` to use a `DIMER` type interaction. This interaction really means that the two binding sites must have a specific orientation relative to each other.


In general, lines like this have the form

```
<name1>	<name2>	<coop_type>	<distance_threshold>	<any>	<other>	<parameters>	<specific>	<to>	<the>	<interaction>
```

Where the number of parameters after the `distance_threshold` depend on the type of interaction being requested.

## Types of special cooperativities

### DIMER interaction

The DIMER interaction specifies a particular orientaiton configuration.

```
<name1>	<name2>	DIMER <distance> <strand_a> <strand_b>
```

requires the two sites to be a name1 site on strand_a and a name2 site on strand_b. The the reverse configuration on the opposite strands (reverse complement sequence) will also work.

HOMODIMERS
If you use the DIMER interaction between a TF and itself, then you must have opposite strands for strand_a and strand_b. If you would like a long chain interaction, just use a normal interaction.

NOTE
Despite its name, the DIMER interaction can be used in other ways, such as to create complexes of several TFs that must be in a particular order.
