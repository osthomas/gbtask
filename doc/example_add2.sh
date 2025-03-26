gbtask modify --pad-left 100 --pad-right AAAAGGGGGCCCCCTTTTT \
    --feature 10 20 1 exon --feature 11 19 1 CDS \
    -q all mygroupingqualifier mygene \
    -i ex1.gb > ex2.gb
cat ex2.gb
