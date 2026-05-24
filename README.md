# Bit-Pep: Multi-Proteome Peptide Classification

Adapting [bit-pop](https://github.com/mladenpop-oss/bit-pop/) to map Peptide(s) to Proteome(s) as [bit-pep](https://github.com/animesh/bit-pep/) , running it is simple as installing cargo/rust `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`, adding to the path, for bash/ubuntu `source $HOME/.bashrc` and running `cargo run --bin bit-pep -- run-prot <proteome in fasta format> -p <peptides list one per line in txt or fasta format> -j <number of CPU to use> -m <RAM in GB>` for example to map peptides from MaxQuant HeLa DIA [results](https://fuzzylife.substack.com/p/proteomics-data-processing-with-maxquant) using 2 CPU (-j 2) and 4GB RAM (-m 2)

```
# install cargo/rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# update path 
source $HOME/.bashrc
# run bit-pep
cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta -p peptides.txt -j 2 -m 4
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.76s
     Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p peptides.txt -j 2 -m 4`
bit-pep RunProt  (FM-index parallel peptide->proteome search)
===============================================================
  Threads : 2 (of 2 logical CPUs)
  Memory  : 4.0 GB  ->  batch size: 33554432 peptides

[1/4] Resolving proteome...
  uniprot_sprot.fasta (274.2 MB)

[2/4] Indexing proteins...
  Loading: [00:00:18 ████████████████████████████████████████ 574627/574627]                                          574627 proteins indexed in 75.0s  (RAM: 3430 MB)

[3/4]+[4/4] Streaming peptides and mapping (2 threads, batch=33554432)...
  done -- 58920 submitted, 54626 mapped (92.7%), 3.1s                                                               
================================================================
  Output  : peptides.pep.tsv
  Mapping : 3.09s  |  Total: 78.12s  |  Peak RAM: 3434 MB

-----------------------------------------
Peptides submitted : 58920
Mapped             : 54626 (92.7%)
  unique            : 25225
  shared (>1 prot.) : 208666
  cross-proteome    : 0
Unmapped           : 4294
-----------------------------------------

Peptide -> Protein -> Species  (top 20 by total peptides)
  --------------------------------------  --------  --------  --------  ------  --------  --------  ------
  Organism                                  Unique    Shared     Total  %total  ProtsHit  TotalProt  %Prots
  --------------------------------------  --------  --------  --------  ------  --------  --------  ------
  Homo sapiens                               25107     50456     75563   128.2      6852     20431    33.5
  Mus musculus                                  49     29936     29985    50.9      4656     17252    27.0
  Rattus norvegicus                             18     16870     16888    28.7      2333      8226    28.4
  Bos taurus                                     4     12063     12067    20.5      1651      6052    27.3
  Pongo abelii                                   7      8660      8667    14.7       898      2218    40.5
  Gallus gallus                                  3      5230      5233     8.9       634      2314    27.4
  Xenopus laevis                                 0      4372      4372     7.4       752      3514    21.4
  Sus scrofa                                     0      3933      3933     6.7       439      1462    30.0
  Macaca fascicularis                            1      3758      3759     6.4       417      1176    35.5
  Oryctolagus cuniculus                          1      3042      3043     5.2       348       979    35.5
  Canis lupus familiaris                         0      2356      2356     4.0       251       857    29.3
  Pan troglodytes                                0      2273      2273     3.9       224       692    32.4
  Arabidopsis thaliana                           1      2264      2265     3.8       533     16418     3.2
  Drosophila melanogaster                        1      2206      2207     3.7       366      3868     9.5
  Danio rerio                                    1      2039      2040     3.5       502      3369    14.9
  Xenopus tropicalis                             0      1659      1659     2.8       326      1713    19.0
  Caenorhabditis elegans                         0      1555      1555     2.6       294      4499     6.5
  Cricetulus griseus                             0      1526      1526     2.6       105       248    42.3
  Dictyostelium discoideum                       0      1159      1159     2.0       196      4163     4.7
  Mesocricetus auratus                           0      1056      1056     1.8       102       277    36.8
  --------------------------------------  --------  --------  --------  ------  --------  --------  ------
  ... and 2663 more organisms

  Unique/Shared/Total = peptide counts; %total = % of submitted
  ProtsHit/TotalProt = proteins hit / proteins in DB; %Prots = coverage
```

notice `4294` peptides are `Unmapped` which is mainly due to missing isofom/peptides[LAQPGFPSGGPGGTR,SPIAAARCR,...] in [uniprot_sprot.fasta](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) and deleted protein/peptides[ELQQVTAGEAASIH,QVCQVPASR,...] in [current human proteome](https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&includeIsoform=true&query=%28%28proteome%3AUP000005640%29%29) leading to `Unmapped           : 4294` verified by [including](https://zenodo.org/records/20344308) human proteome from [2024](https://github.com/user-attachments/files/16439166/mqpar.xml.txt)
```
wget "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&includeIsoform=true&query=%28%28proteome%3AUP000005640%29%29" -O human.fasta
wget https://zenodo.org/records/20344308/files/uniprotkb_proteome_UP000005640_2024_04_18.fasta 
cat uniprotkb_proteome_UP000005640_2024_04_18.fasta >> human.fasta
wget "https://zenodo.org/records/14557756/files/proteinGroups.txt"
grep -vE "Peptide|CON__|REV__" proteinGroups.txt| awk -F '\t' '{print $92}'  | sed 's/;/\n/g'  > peptides.txt 
cargo run --bin bit-pep -- run-prot human.fasta  -p peptides.txt -j 12 
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 3.56s
        Running `target/debug/bit-pep run-prot human.fasta -p peptides.txt -j 12`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================
    Threads : 12 (of 12 logical CPUs)
    Memory  : 67.5 GB  ->  batch size: 565961632 peptides

    [1/4] Resolving proteome...
    human.fasta (106.3 MB)

    [2/4] Indexing proteins...
    Loading: [00:00:07 ████████████████████████████████████████ 210230/210230]                                                         210230 proteins indexed in 49.2s  (RAM: 1392 MB)

    [3/4]+[4/4] Streaming peptides and mapping (12 threads, batch=565961632)...
    done -- 58920 submitted, 58920 mapped (100.0%), 6.7s                                                                             
    ================================================================
    Output  : peptides.pep.tsv
    Mapping : 6.70s  |  Total: 55.85s  |  Peak RAM: 1401 MB

    -----------------------------------------
    Peptides submitted : 58920
    Mapped             : 58920 (100.0%)
    unique            : 23403
    shared (>1 prot.) : 657336
    cross-proteome    : 0
    Unmapped           : 0
    -----------------------------------------

    Peptide -> Protein -> Species  (top 20 by total peptides)
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Organism                                  Unique    Shared     Total  %total  ProtsHit  TotalProt  %Prots
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Homo sapiens                               23403    657336    680739  1155.4     35721    210230    17.0
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------

    Unique/Shared/Total = peptide counts; %total = % of submitted
    ProtsHit/TotalProt = proteins hit / proteins in DB; %Prots = coverage
```

### position of the hit is written as tab-separated text file in `<proteome>.pep.tsv` 

for example position of hits for [2 selected peptides](https://raw.githubusercontent.com/animesh/bit-pep/refs/heads/master/pep2Q9UPZ3.txt) from [Q9UPZ3](https://www.uniprot.org/uniprotkb/Q9UPZ3/entry) [sequences](https://rest.uniprot.org/uniprotkb/Q9UPZ3.fasta) with an old PERL [code](https://raw.githubusercontent.com/animesh/scripts/079de9ce5a6d408d7508fb276f6ad248c58b99d7/pep2protmap.pl) 
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget https://raw.githubusercontent.com/animesh/scripts/079de9ce5a6d408d7508fb276f6ad248c58b99d7/pep2protmap.pl
perl pep2protmap.pl uniprot_sprot.fasta pep2Q9UPZ3.txt 1

    Read# 574627 sequences from uniprot_sprot.fasta

    Opening peptide list from pep2Q9UPZ3.txt

    >sp|Q9UPZ3|HPS5_HUMAN   569;
    >sp|Q9UPZ3|HPS5_HUMAN   111;

    Processed 2 Sequences
    Found 2 Matches

 cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta  -p pep2Q9UPZ3.txt
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.65s
     Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p pep2Q9UPZ3.txt`
bit-pep RunProt  (FM-index parallel peptide->proteome search)
===============================================================
  Threads : 1 (of 2 logical CPUs)
  Memory  : 1.0 GB  ->  batch size: 8388608 peptides

[1/4] Resolving proteome...
  uniprot_sprot.fasta (274.2 MB)

[2/4] Indexing proteins...
  Loading: [00:00:18 ████████████████████████████████████████ 574627/574627]                                          574627 proteins indexed in 75.0s  (RAM: 3429 MB)

[3/4]+[4/4] Streaming peptides and mapping (1 threads, batch=8388608)...
  done -- 2 submitted, 2 mapped (100.0%), 0.5s                                                                      
================================================================
  Output  : pep2Q9UPZ3.pep.tsv
  Mapping : 0.46s  |  Total: 75.51s  |  Peak RAM: 3429 MB

-----------------------------------------
Peptides submitted : 2
Mapped             : 2 (100.0%)
  unique            : 2
  shared (>1 prot.) : 0
  cross-proteome    : 0
Unmapped           : 0
-----------------------------------------

Peptide -> Protein -> Species  (top 20 by total peptides)
  --------------------------------------  --------  --------  --------  ------  --------  --------  ------
  Organism                                  Unique    Shared     Total  %total  ProtsHit  TotalProt  %Prots
  --------------------------------------  --------  --------  --------  ------  --------  --------  ------
  Homo sapiens                                   2         0         2   100.0         1     20431     0.0
  --------------------------------------  --------  --------  --------  ------  --------  --------  ------

  Unique/Shared/Total = peptide counts; %total = % of submitted
  ProtsHit/TotalProt = proteins hit / proteins in DB; %Prots = coverage

cat pep2Q9UPZ3.pep.tsv 
    peptide_id      sequence        protein_acc     protein_name    proteome_id     start   end     score   mismatches status
    p0      PELRGDEQSCEEDVSSDTCPK   Q9UPZ3  HPS5_HUMAN BLOC-2 complex member HPS5 [HPS5] Homo sapiens       uniprot_sprot.fasta 570     590     1.0000  0       unique
    p1      PEQMYVSSEHK     Q9UPZ3  HPS5_HUMAN BLOC-2 complex member HPS5 [HPS5] Homo sapiens       uniprot_sprot.fasta112      122     1.0000  0       unique
```

### for large peptide list ~100M

increase memory up to to 24GB (-m 24) otherwise CRASH is iminent, like for [all [6-7](https://en.wikipedia.org/wiki/6-7) lengths] over uniprot to fully check for overlapping peptides by [Chopping Proteins to Peptides](https://fuzzylife.substack.com/p/chopping-proteins-to-peptides), notice the peak mem usage is ~12GB though, maybe it is 6-7 the issue?
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget https://raw.githubusercontent.com/animesh/scripts/d08e6ffb9c9fdc23f37fd12cdd280966d4c9e46e/pepCleave.pl
perl pepCleave.pl uniprot_sprot.fasta 6 7
grep -v "^>" uniprot_sprot.fasta.len6to7.fasta > uniprot_sprot.fasta.len6to7.txt
wc uniprot_sprot.fasta.len6to7.*
  246986978  2411497827 23586051901 uniprot_sprot.fasta.len6to7.fasta
  123493489   123493489   953910896 uniprot_sprot.fasta.len6to7.txt
  370480467  2534991316 24539962797 total
cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta -p uniprot_sprot.fasta.len6to7.txt -j 12 -m 48
    Compiling bit-pep v0.2.0 (/mnt/z/Download/bit-pep)
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 19.61s
        Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p uniprot_sprot.fasta.len6to7.txt -j 12 -m 48`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================
    Threads : 12 (of 12 logical CPUs)
    Memory  : 48.0 GB  ->  batch size: 402653184 peptides

    [1/4] Resolving proteome...
    uniprot_sprot.fasta (274.2 MB)

    [2/4] Indexing proteins...
    Loading: [00:00:19 ████████████████████████████████████████ 574627/574627]                                                         574627 proteins indexed in 117.0s  (RAM: 3429 MB)

    [3/4]+[4/4] Streaming peptides and mapping (12 threads, batch=402653184)...
    ⠁ 121000000 searched, 120999999 mapped, RAM 16670 MB, 5886.8s                                                                      
    done -- 123493489 submitted, 123493489 mapped (100.0%), 6004.8s                                                                  
    ================================================================
    Output  : uniprot_sprot.fasta.len6to7.pep.tsv
    Mapping : 6004.84s  |  Total: 6121.85s  |  Peak RAM: 11018 MB

    -----------------------------------------
    Peptides submitted : 123493489
    Mapped             : 123493489 (100.0%)
    unique            : 68518651
    shared (>1 prot.) : 54974838
    cross-proteome    : 0
    Unmapped           : 0
    -----------------------------------------

    Peptide -> Protein -> Species  (top 20 by total peptides)
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Organism                                  Unique    Shared     Total  %total  ProtsHit  TotalProt  %Prots
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Homo sapiens                             3972913  18633102  22606015    18.3     20427     20431   100.0
    Mus musculus                             3075452  16358493  19433945    15.7     17251     17252   100.0
    Arabidopsis thaliana                     4902809   9793502  14696311    11.9     16413     16418   100.0
    Rattus norvegicus                         780436   7595196   8375632     6.8      8226      8226   100.0
    Saccharomyces cerevisiae (strain AT...   2114681   3863170   5977851     4.8      6733      6733   100.0
    Drosophila melanogaster                  1788004   2990340   4778344     3.9      3868      3868   100.0
    Bos taurus                                639628   4135212   4774840     3.9      6052      6052   100.0
    Caenorhabditis elegans                   1910678   2801063   4711741     3.8      4499      4499   100.0
    Schizosaccharomyces pombe (strain 9...   1943759   2766388   4710147     3.8      5129      5129   100.0
    Dictyostelium discoideum                 1700704   2891397   4592101     3.7      4163      4163   100.0
    Oryza sativa subsp. japonica              909290   2937931   3847221     3.1      4197      4197   100.0
    Danio rerio                              1054844   2209411   3264255     2.6      3369      3369   100.0
    Xenopus laevis                            670087   2530713   3200800     2.6      3514      3514   100.0
    Escherichia coli (strain K12)             498285   2227726   2726011     2.2      4531      4531   100.0
    Bacillus subtilis (strain 168)            874821   1540350   2415171     2.0      4191      4191   100.0
    Gallus gallus                             513067   1742942   2256009     1.8      2313      2314   100.0
    Pongo abelii                               54298   1897918   1952216     1.6      2218      2218   100.0
    Mycobacterium tuberculosis (strain ...     97687   1497889   1595576     1.3      2338      2338   100.0
    Xenopus tropicalis                        281088   1204105   1485193     1.2      1713      1713   100.0
    Mycobacterium tuberculosis (strain ...      1407   1285584   1286991     1.0      1899      1899   100.0
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    ... and 14821 more organisms

    Unique/Shared/Total = peptide counts; %total = % of submitted
    ProtsHit/TotalProt = proteins hit / proteins in DB; %Prots = coverage
```


check with random peptides generated per protein from 10% of uniprot peptides ranging from 5 to 35 with bit-pep default using half of logical CPU available
```
python generate_random_peptides.py 
    Sampled 57455/574551 proteins (10.0%)
    Generated 57455 peptides from sampled proteins in uniprot_sprot.fasta
    Wrote peptide FASTA: uniprot_sprot_0.1_5_35_42.fasta
    Running: cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta -p uniprot_sprot_0.1_5_35_42.fasta
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 3.78s
        Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p uniprot_sprot_0.1_5_35_42.fasta`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================
    Threads : 6 (of 12 logical CPUs)
    Memory  : 67.1 GB  ->  batch size: 562943328 peptides

    [1/4] Resolving proteome...
    uniprot_sprot.fasta (274.2 MB)

    [2/4] Indexing proteins...
    Loading: [00:00:19 ████████████████████████████████████████ 574627/574627]                                                         574627 proteins indexed in 117.1s  (RAM: 3429 MB)

    [3/4]+[4/4] Streaming peptides and mapping (6 threads, batch=562943328)...
    done -- 57455 submitted, 57455 mapped (100.0%), 10.5s                                                                            
    ================================================================
    Output  : uniprot_sprot_0.1_5_35_42.pep.tsv
    Mapping : 10.49s  |  Total: 127.56s  |  Peak RAM: 3465 MB

    -----------------------------------------
    Peptides submitted : 57455
    Mapped             : 57455 (100.0%)
    unique            : 27748
    shared (>1 prot.) : 899722
    cross-proteome    : 0
    Unmapped           : 0
    -----------------------------------------

    Peptide -> Protein -> Species  (top 20 by total peptides)
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Organism                                  Unique    Shared     Total  %total  ProtsHit  TotalProt  %Prots
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Homo sapiens                                1100     27065     28165    49.0     12986     20431    63.6
    Mus musculus                                 925     22524     23449    40.8     11099     17252    64.3
    Arabidopsis thaliana                        1355     14809     16164    28.1      9277     16418    56.5
    Rattus norvegicus                            346     10164     10510    18.3      5274      8226    64.1
    Bos taurus                                   312      6426      6738    11.7      3658      6052    60.4
    Escherichia coli (strain K12)                210      6217      6427    11.2      2700      4531    59.6
    Saccharomyces cerevisiae (strain AT...       543      5772      6315    11.0      3554      6733    52.8
    Oryza sativa subsp. japonica                 273      5747      6020    10.5      2866      4197    68.3
    Drosophila melanogaster                      314      5682      5996    10.4      2458      3868    63.5
    Caenorhabditis elegans                       409      4298      4707     8.2      2485      4499    55.2
    Schizosaccharomyces pombe (strain 9...       467      4191      4658     8.1      2731      5129    53.2
    Dictyostelium discoideum                     363      4168      4531     7.9      2233      4163    53.6
    Escherichia coli O157:H7                      12      4511      4523     7.9      1489      2047    72.7
    Escherichia coli O6:H1 (strain CFT0...         6      4137      4143     7.2      1302      1702    76.5
    Shigella flexneri                             10      3995      4005     7.0      1278      1696    75.4
    Salmonella typhimurium (strain LT2 ...        34      3846      3880     6.8      1336      1831    73.0
    Xenopus laevis                               219      3656      3875     6.7      2102      3514    59.8
    Danio rerio                                  285      3282      3567     6.2      1964      3369    58.3
    Salmonella typhi                               6      3400      3406     5.9      1077      1349    79.8
    Shigella sonnei (strain Ss046)                 1      3348      3349     5.8       824       879    93.7
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    ... and 12424 more organisms

    Unique/Shared/Total = peptide counts; %total = % of submitted
    ProtsHit/TotalProt = proteins hit / proteins in DB; %Prots = coverage
    Wrote metadata CSV: uniprot_sprot_0.1_5_35_42.metadata.csv
    Wrote verified report: uniprot_sprot_0.1_5_35_42.verified.txt

cat uniprot_sprot_0.1_5_35_42.verified.txt
    Total peptides: 57455
    Verified peptides: 57455
    Unique peptides: 27548
    Shared peptides: 29907
    Missing/mismatched peptides: 0

awk -F ',' '{print $NF}' uniprot_sprot_0.1_5_35_42.metadata.csv | sort | uniq -c
  29907 shared
  27548 unique
      1 verified
```

### test sweep random seeds from 1 to 100 with peptides ranging from 5-35 per protein from ~1/3 (0.33) of uniprot and pops-up [Dictyostelium discoideum](https://en.wikipedia.org/wiki/Dictyostelium_discoideum) seems to possesses a massive number of highly conserved eukaryotic proteins!

```
for seed in $(seq 1 100); do echo $seed;
  python3 generate_random_peptides.py \
    --fasta uniprot_sprot.fasta \
    --min-len 5 \
    --max-len 35 \
    --sample-fraction 0.33 \
    --seed "$seed"
done
    #seed 1 
    #... to 99 not shown
    100
    Sampled 189601/574551 proteins (33.0%)
    Generated 189601 peptides from sampled proteins in uniprot_sprot.fasta
    Wrote peptide FASTA: uniprot_sprot_0.33_5_35_100.fasta
    Running: cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta -p uniprot_sprot_0.33_5_35_100.fasta
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 3.02s
        Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p uniprot_sprot_0.33_5_35_100.fasta`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================
    Threads : 6 (of 12 logical CPUs)
    Memory  : 67.2 GB  ->  batch size: 563494240 peptides

    [1/4] Resolving proteome...
    uniprot_sprot.fasta (274.2 MB)

    [2/4] Indexing proteins...
    Loading: [00:00:19 ████████████████████████████████████████ 574627/574627]                                                         574627 proteins indexed in 118.7s  (RAM: 3429 MB)

    [3/4]+[4/4] Streaming peptides and mapping (6 threads, batch=563494240)...
    done -- 189601 submitted, 189601 mapped (100.0%), 33.5s                                                                          
    ================================================================
    Output  : uniprot_sprot_0.33_5_35_100.pep.tsv
    Mapping : 33.46s  |  Total: 152.14s  |  Peak RAM: 3533 MB

    -----------------------------------------
    Peptides submitted : 189601
    Mapped             : 189601 (100.0%)
    unique            : 89888
    shared (>1 prot.) : 3185051
    cross-proteome    : 0
    Unmapped           : 0
    -----------------------------------------

    Peptide -> Protein -> Species  (top 20 by total peptides)
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Organism                                  Unique    Shared     Total  %total  ProtsHit  TotalProt  %Prots
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    Dictyostelium discoideum                    1231    130917    132148    69.7      3727      4163    89.5
    Homo sapiens                                3694     97499    101193    53.4     18988     20431    92.9
    Mus musculus                                3018     83056     86074    45.4     16143     17252    93.6
    Arabidopsis thaliana                        4503     53827     58330    30.8     14778     16418    90.0
    Rattus norvegicus                           1047     37303     38350    20.2      7704      8226    93.7
    Drosophila melanogaster                     1053     26015     27068    14.3      3574      3868    92.4
    Saccharomyces cerevisiae (strain AT...      1810     22827     24637    13.0      5821      6733    86.5
    Bos taurus                                   967     22578     23545    12.4      5563      6052    91.9
    Oryza sativa subsp. japonica                 879     22351     23230    12.3      4003      4197    95.4
    Escherichia coli (strain K12)                608     21247     21855    11.5      4046      4531    89.3
    Caenorhabditis elegans                      1313     15654     16967     8.9      4038      4499    89.8
    Schizosaccharomyces pombe (strain 9...      1485     14310     15795     8.3      4452      5129    86.8
    Escherichia coli O157:H7                      29     15411     15440     8.1      1969      2047    96.2
    Xenopus laevis                               718     13496     14214     7.5      3265      3514    92.9
    Escherichia coli O6:H1 (strain CFT0...        19     14136     14155     7.5      1645      1702    96.7
    Shigella flexneri                             30     13636     13666     7.2      1634      1696    96.3
    Danio rerio                                  922     12370     13292     7.0      3116      3369    92.5
    Salmonella typhimurium (strain LT2 ...       120     12802     12922     6.8      1755      1831    95.8
    Escherichia coli O9:H4 (strain HS)             4     11445     11449     6.0       847       847   100.0
    Escherichia coli O139:H28 (strain E...         8     11372     11380     6.0       838       838   100.0
    --------------------------------------  --------  --------  --------  ------  --------  --------  ------
    ... and 14225 more organisms

    Unique/Shared/Total = peptide counts; %total = % of submitted
    ProtsHit/TotalProt = proteins hit / proteins in DB; %Prots = coverage
    Wrote metadata CSV: uniprot_sprot_0.33_5_35_100.metadata.csv
    Wrote verified report: uniprot_sprot_0.33_5_35_100.verified.txt

cat uniprot_sprot_0.33*veri* | grep ssing | sort | uniq -c
    100 Missing/mismatched peptides: 0

awk -F ',' '{print $NF}' uniprot_sprot_0.33*.metadata.csv | sort | uniq -c
9974531 shared
8985569 unique
    100 verified
```

# change summary with @copilot

## 1. AA encoding

### `src/lib.rs`
```rust
pub fn encode_sequence_aa(seq: &str) -> Vec<u8> {
    seq.bytes().filter_map(|c| {
        let code = crate::amino::encode_aa(c);
        if code >= 1 { Some(code) } else { None }
    }).collect()
}
```

### `src/amino.rs`
```rust
pub const ALPHA_SIZE: usize = 27;

pub fn encode_aa(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' =>  1,
        b'C' =>  2,
        b'D' =>  3,
        b'E' =>  4,
        b'F' =>  5,
        b'G' =>  6,
        b'H' =>  7,
        b'I' =>  8,
        b'K' =>  9,
        b'L' => 10,
        b'M' => 11,
        b'N' => 12,
        b'P' => 13,
        b'Q' => 14,
        b'R' => 15,
        b'S' => 16,
        b'T' => 17,
        b'V' => 18,
        b'W' => 19,
        b'Y' => 20,
        b'U' => 21,
        b'O' => 22,
        b'B' => 23,
        b'J' => 24,
        b'Z' => 25,
        b'X' => 26,
        _    => 26,
    }
}

pub fn decode_aa(code: u8) -> u8 {
    match code {
        1..=20 => AA_ORDER[(code - 1) as usize],
        21 => b'U',
        22 => b'O',
        23 => b'B',
        24 => b'J',
        25 => b'Z',
        _  => b'X',
    }
}

pub fn encode_sequence(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        if c == b'*' { break; }
        out.push(encode_aa(c));
    }
    out
}
```

## 2. Protein mode flag and accessors

### `src/lib.rs`
```rust
protein_mode: bool,  // true = AA encoding, no reverse complement
```

```rust
pub fn set_protein_mode(&mut self, v: bool) { self.protein_mode = v; }
pub fn is_protein_mode(&self) -> bool { self.protein_mode }
```

## 3. Protein genome loading

### `src/lib.rs`
```rust
pub fn add_genome(&mut self, name: &str, sequence: &str) -> u32 {
    let genome_id = self.genomes.len() as u32;
    let encoded = if self.protein_mode { encode_sequence_aa(sequence) } else { encode_sequence(sequence) };
    self.genome_names.insert(genome_id, name.to_string());
    self.genomes.insert(genome_id, encoded.clone());
    genome_id
}
```

## 4. Read encoding chooses AA in protein mode

### `src/lib.rs`
```rust
fn encode_read(&self, read: &str) -> Vec<u8> {
    if self.protein_mode { encode_sequence_aa(read) } else { encode_sequence(read) }
}
```

## 5. No reverse complement for proteins

### `src/lib.rs`
```rust
pub fn map_read_with_mode(
    &self,
    read: &str,
    mode: AlignMode,
    context_window: usize,
) -> Vec<MappingResult> {
    let forward_results = self.map_read_orientation(read, mode, context_window, false);
    if self.protein_mode { return forward_results; }  // proteins have no strand
    let rc_read = reverse_complement(read);
    let rc_results = self.map_read_orientation(&rc_read, mode, context_window, true);
    ...
}
```

## 6. Protein-specific scoring and alignment branch

### `src/lib.rs`
```rust
if self.protein_mode {
    let start = position as isize - anchor_read_offset as isize;
    if start < 0 { continue; }
    let start = start as usize;
    let end = start + read_len;
    if end > genome.len() { continue; }
    let region = &genome[start..end];
    if self.fuzzy_mismatches == 0 || matches!(self.fuzzy_method, FuzzyMethod::None) {
        if region == encoded.as_slice() {
            let cigar = format!("{}M", read_len);
            scored.push((genome_id, start as u64, 1.0f64, cigar));
        }
    } else {
        let mismatches = region.iter().zip(encoded.iter())
            .filter(|(a, b)| a != b)
            .count();
        if mismatches <= self.fuzzy_mismatches {
            let score = (read_len - mismatches) as f64 / read_len as f64;
            let cigar = region.iter().zip(encoded.iter())
                .map(|(a, b)| if a == b { 'M' } else { 'X' })
                .collect::<String>();
            scored.push((genome_id, start as u64, score, cigar));
        }
    }
    continue;
}
```

## 7. Protein mode uses all positions and disables subsampling

### `src/lib.rs`
```rust
let pos_cap = if self.protein_mode { usize::MAX } else { 500 };
```

```rust
let positions: Vec<(u32, u64)> = if self.protein_mode {
    raw_positions
} else if raw_positions.len() > 100 {
    let stride = raw_positions.len() / 100;
    raw_positions.into_iter().step_by(stride).collect()
} else {
    raw_positions
};
```

## 8. Protein mode disables k-mer repeat filtering

### `src/lib.rs`
```rust
let repeat_thresh = if self.protein_mode { usize::MAX } else { DEFAULT_REPEAT_THRESHOLD };
```

## 9. Do not truncate protein results

### `src/lib.rs`
```rust
if !self.protein_mode {
    scored.truncate(50);
}
```

## 10. Peptide-specific APIs

### `src/lib.rs`
```rust
pub fn find_peptide_exact(&self, peptide: &str) -> Vec<(u32, u64)> {
    let fm = match &self.fm_index {
        Some(f) => f,
        None    => return Vec::new(),
    };
    let encoded = encode_sequence_aa(peptide);
    fm.find_positions(&encoded, usize::MAX)
}
```

```rust
pub fn find_peptide_fuzzy(&self, peptide: &str, max_mismatches: usize) -> Vec<(u32, u64, u32)> {
    let fm = match &self.fm_index {
        Some(f) => f,
        None    => return Vec::new(),
    };
    let encoded = encode_sequence_aa(peptide);
    let pep_len = encoded.len();
    if pep_len == 0 { return Vec::new(); }

    let n_seeds  = max_mismatches + 1;
    let seed_len = (pep_len / n_seeds).max(1);

    let mut candidates: std::collections::HashSet<(u32, u64)> = std::collections::HashSet::new();

    for s in 0..n_seeds {
        let seed_start = s * seed_len;
        let seed_end   = ((s + 1) * seed_len).min(pep_len);
        let seed       = &encoded[seed_start..seed_end];
        for (genome_id, seed_pos) in fm.find_positions(seed, usize::MAX) {
            if seed_pos < seed_start as u64 { continue; }
            let cand_start = seed_pos - seed_start as u64;
            candidates.insert((genome_id, cand_start));
        }
    }

    let mut results = Vec::new();
    for (genome_id, pos0) in candidates {
        if let Some(genome) = self.genomes.get(&genome_id) {
            let end = pos0 as usize + pep_len;
            if end > genome.len() { continue; }
            let region = &genome[pos0 as usize..end];
            let mismatches = region.iter().zip(encoded.iter())
                .filter(|(a, b)| a != b).count() as u32;
            if mismatches <= max_mismatches as u32 {
                results.push((genome_id, pos0, mismatches));
            }
        }
    }
    results
}
```

## 11. Protein runtime entry point in `src/main.rs`

```rust
let mut bp = bit_pop::BitPop::new(5);   // k=5 optimal for protein 5-mers
bp.set_protein_mode(true);
```

## 12. UniProt/proteome support in `src/main.rs`

`src/main.rs` includes protein-specific runtime helpers:
- `parse_uniprot_header(header: &str) -> (String, String, String, String)` parses UniProt FASTA headers into accession, protein name, organism, and gene.
- `resolve_prot_fasta(proteome: &str, force: bool) -> Result<PathBuf, String>` resolves `sprot`, `trembl`, `isoforms`, local FASTA paths, or cached UniProt proteome IDs.

These helpers are used by `cmd_run_prot` to build the protein database and metadata for TSV output.

## 13. Protein/peptide integration tests

### `tests/peptide_map_integration.rs`
- downloads a UniProt proteome via `bit_pop::uniprot::fetch_proteome`
- parses FASTA headers and sequences
- encodes proteins with `bit_pop::amino::encode_sequence`
- digests proteins with trypsin-style cleavage
- parses peptide input with `bit_pop::peptide_input::parse_peptide_str`
- exercises peptide uniqueness and valid AA encoding

## 14. Supporting protein-specific modules

- `src/peptide_input.rs` — peptide file parsing for protein mode
- `src/pep_output.rs` — peptide/protein TSV hit reporting and summary status
- `src/uniprot.rs` — UniProt proteome download/search support and proteome entry display
- `src/amino.rs` — amino acid alphabet encoding/decoding and packing

## 13. Protein runtime entrypoint in `src/main.rs`

The protein mode entry point is the `RunProt` subcommand.
It:
- defines `RunProtArgs` for `--peptides`, `--fuzzy-mismatches`, `-j`, `-m`, `--force`
- resolves a proteome path or UniProt UPID via `resolve_prot_fasta()`
- builds `BitPop` with `BitPop::new(5)` and calls `bp.set_protein_mode(true)`
- reads proteins from FASTA and stores metadata (accession, name, organism, gene)
- streams peptides in batches from the input file
- uses `bp.find_peptide_exact()` or `bp.find_peptide_fuzzy()` per peptide
- assembles `PepHit` rows and writes TSV output

## 14. Peptide input parser in `src/peptide_input.rs`

This module supports:
- plain mode: one peptide per line or space-separated peptides
- FASTA mode: sequences under `>` headers
- validation for ASCII amino-acid characters and length
- stop-codon trimming (`*`)
- warnings for invalid or too-long sequences

It exposes:
- `Peptide { id, sequence }`
- `parse_peptide_file(path)`
- `parse_peptide_str(input)`

## 15. Peptide TSV output in `src/pep_output.rs`

This module contains:
- `PepHit` record type with fields `peptide_id`, `sequence`, `protein_acc`, `protein_name`, `proteome_id`, `start`, `end`, `score`, `mismatches`, `status`
- `HitStatus` enum: `Unique`, `Shared`, `CrossProteome`, `Unmapped`
- `assign_status(hits)` to classify peptide hits
- `write_tsv(hits, path)` to write output
- `print_summary(hits, n_peptides)` for a summary report

## Additional protein-specific references exist in non-runtime file:
- `src/MIGRATION.md` — migration notes and source-file changes


# TBD 

FM-index gives compaction, exact-matching and scalable search and proteins have high alphabet entropy so 5-mer already has strong selectivity: 26^5 \approx 1.19 \times 10^7 so unlike DNA, candidate explosion is naturally limited but but there are many opportunities for improvemements for sure, like 
* index structure tuning
* SIMD verification
* seed strategy
* memory layout
* cache efficiency
* ...?
## FM-index maybe not beat hash indexing?
For short exact peptide lookup, FM-index is elegant but may not be fastest because proteins have high entropy, a direct k-mer hash table may outperform FM traversal. For example k=6:26^6 \approx 3.09 \times 10^8 but in practice only a tiny fraction exist in real proteomes. So sparse hash maps become feasible with lookup advantages like :
* O(1) average lookup
* very cache friendly
* simpler mismatch neighborhood expansion
* easier SIMD integration
For modern servers with abundant RAM, hash-based indexes often win in throughput. Need to benchmark!

## mismatch verification should become SIMD? 
```rust id="th0u79"
zip()
filter()
count()
```
For millions of peptides this becomes dominant as we saw less than -m 24 leads to crash and AA we already encode residues into compact integers which is ideal for SIMD. Potential AVX2 approach:
* load 32 bytes
* compare vectors
* movemask
* popcount mismatches
* early terminate after >2 mismatches
... could produce very large speedups? 

## k=5 may not actually be optimal
currently use:
```rust id="jlwmv9"
BitPop::new(5)
```
allowing 1–2 mismatches, fixed k=5 may actually be suboptimal.
need to empirically benchmark:
* k=5
* k=6
* k=7
* ... spaced seeds?

## seed-and-extend strategy can improve?
Current fuzzy search:
```rust id="jlwmby"
n_seeds = mismatches + 1
```
maybe better to:
* spaced seeds
* entropy-aware seeds
* rare-seed prioritization
* ...?

## deduplication may become bottleneck
```rust id="vn1y6m"
HashSet<(u32, u64)>
```
is potentially expensive. For high-throughput peptide search:
* hashing overhead
* allocator churn
* poor cache locality
can dominate sorted vectors + dedup may outperform hashsets substantially, needs benchmarking for sure!

## ambiguous residue handling is probably too permissive
unknown residues map to X so all unknowns treated equivalently and this may inflate false positives.  Maybe:
* configurable ambiguity policy
* penalized ambiguity matching
* reject high-X queries
* ...?

## not to mention the biggest missing feature is PTM-aware matching!!


## Availability

Source code available under the MIT License.

## Citation

```bibtex
@software{popovic_2026_bitpop,
  author = {Popovi{\'c}, Mladen},
  title = {Bit-Pop: A Proof-of-Concept Tool for Multi-Genome DNA Read Classification},
  year = {2026},
  doi = {10.5281/zenodo.20043593},
  url = {https://github.com/mladenpop-oss/bit-pop}
}
```

Or plain text:

> Popović, M. (2026). Bit-Pop: A Proof-of-Concept Tool for Multi-Genome DNA Read Classification. https://doi.org/10.5281/zenodo.20043593

## License

MIT License
