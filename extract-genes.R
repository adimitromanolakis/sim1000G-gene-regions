library(sim1000G)

## Script to extract genes from 1000 genomes 20130502 data

## 1000 genomes data were downloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

## Genes file was extracted from ensembl version GRCH37:
##   https://grch37.ensembl.org/biomart/martview


# To run in parallel do something like:
#  for i in 1 2 3 4 5 6; do Rscript extract-genes.R $i > log-$i.txt 2>&1 & done
#  for i in 7 8 9 10 11 12 13; do Rscript extract-genes.R $i > log-$i.txt 2>&1 & done
#  for i in 14 15 16 17 18 19 20 21 22; do Rscript extract-genes.R $i > log-$i.txt 2>&1 & done
#


data_dir = "/media/apo/zzz/1000genomes_data/20130502/"
data_dir = "/home/apo/2TB/apostolos/galen-cluster-backup/1000genomes/"

genes_file = file.path(data_dir, "ensembl-genes-grch37.csv")
ped_file   = file.path(data_dir, "integrated_call_samples_v2.20130502.ALL.ped")



stopifnot( file.exists(genes_file) )
stopifnot( file.exists(ped_file) )

if(!exists("genes")) {
    genes <<- read.csv(genes_file, as=T)
    genes <<- genes[ genes$Transcript.type == "protein_coding" , ]
    genes <<- genes[ genes$Chromosome.scaffold.name %in% as.character(1:22) , ]
    table(genes$Chromosome.scaffold.name)
}



extractRegion = function(chrom,
                         bp1,
                         bp2,
                         temp_dir = "/tmp",
                         outfile = "out.vcf",
                         populations = c("CEU","TSI","GBR"),
                         filter = "AC>=1" )
{


    if(!exists("metadata1000genomes")) {
        ped = read.table(ped_file, h=T, sep="\t", as=T)
        metadata1000genomes <<- list(ped = ped)
    }

    ped = metadata1000genomes[["ped"]]

    print(populations)
    selected_indiv = ped$Individual.ID [ ped$Population %in% populations ]

    indiv_file = file.path(temp_dir, "indiv.txt")
    outfile = file.path(temp_dir, outfile)

    cat(selected_indiv, file = indiv_file , sep="\n" )

    region = sprintf("%s:%d-%d", chrom, bp1,bp2)

    vcf = file.path(data_dir, "ALL.chr#CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
    vcf = sub("#CHR", chrom, vcf)

    cmd = sprintf( "bcftools view -S %s --force-samples -r %s %s 2>/tmp/bcftools_log ",
                   indiv_file,
                   region,
                   vcf)

    cmd = sprintf("%s | bcftools filter -i \"%s\" | gzip ", cmd, filter)
    cmd = sprintf("%s > %s ", cmd, outfile)

    cat("--------- Command to run : " , cmd,"\n")
    system(cmd)
}


getGeneCoords = function(gene_name) {



    g = genes[ genes$Gene.name == gene_name, ]
    #print(g)

    chrom = g$Chromosome.scaffold.name[1]
    bp1 = g$Gene.start..bp.[1]
    bp2 = g$Gene.end..bp.[1]

    list(chrom = chrom, bp1 = bp1, bp2 = bp2)
}


extractGeneVcf = function(gene_name,
                          expand_bp = 1000,
                          out_dir = "/tmp",
                          populations = c("CEU","TSI","GBR")) {

    g = getGeneCoords(gene_name)
    cat("Gene coords:::", gene_name, g$chrom, g$bp1, g$bp2, "\n")

    outfile = sprintf("genes-chr%s-%s.vcf.gz", g$chrom, gene_name)

    extractRegion(
        g$chrom,
        g$bp1-expand_bp,
        g$bp2+expand_bp,
        temp_dir = out_dir,
        outfile = outfile,
        filter = 'AC>=1',
        populations = populations
    )
    return(outfile)

}


readGeneVcf = function(gene_name, expand_bp = 1000, out_dir = "/tmp") {

    outfile = extractGeneVcf(gene_name, expand_bp, out_dir)

    vcf = readVCF(
        file.path(out_dir, outfile) ,
        min_maf = 0.01,
        max_maf = 1,
        maxNumberOfVariants = 500
    )

    vcf
}




chr = commandArgs(T)[1]
populations = c("CEU","TSI","GBR")


cat("CHR=",chr,"\n")

{
    genes_chr = genes[genes$Chromosome.scaffold.name == chr,]

    genes_to_extract = unique(genes_chr$Gene.name)
    out_dir = sprintf("/tmp/genes/chr%s",chr)

    system(sprintf("mkdir -p %s", out_dir))

    for(i in 1:length(genes_to_extract)) {
        cat("Extractring gene ",i,"/",length(genes_to_extract),"\n")
    extractGeneVcf(gene_name = genes_to_extract[i], out_dir = out_dir, populations = populations)
    }

}


#extractGeneVcf("FAM72C")
# FAM72C is it long???

getGeneCoords("CKS1B")

exit(1)


geneVcfGitHubLocation = function(gene) {

    #gene = "CKS1B"
    s = genes[genes$Gene.name == gene,]

    chrom = s$Chromosome.scaffold.name[1]

    f="https://raw.githubusercontent.com/adimitromanolakis/sim1000G-gene-regions/main/ceu-tsi-gbr/chr%s/genes-chr%s-%s.vcf.gz"

    f = sprintf(f,chrom,chrom,gene)

    print(f)
    f
}

f="https://raw.githubusercontent.com/adimitromanolakis/sim1000G-gene-regions/main/ceu-tsi-gbr/chr10/genes-chr10-ABCC2.vcf.gz"

vcf = readVCF(geneVcfGitHubLocation("BRCA2"),min_maf = 0)
ped = metadata1000genomes$ped

cat("The vcf file contains populations: ")
table( ped$Population[base::match(vcf$individual_ids,ped$Individual.ID)] )



ped[1,]



vcf$individual_ids







getGeneCoords("KLK9")
genes$L = genes$Gene.end..bp. - genes$Gene.start..bp.
gene_name = "TMEM156"
g = genes[ genes$Gene.name == gene_name, ]



vcf1 = readGeneVcf("KLK10")
vcf2 = readGeneVcf("TMEM156")
vcf3 = readGeneVcf("CNR2")

startMultipleRegionSimulation( list(vcf1,vcf2,vcf3))




ids = generateUnrelatedIndividuals(20)
ids


retrieveGenotypes(52)

fam = newNuclearFamily(1)

x = retrieveGenotypes(fam$gtindex[1])

newFamilyWithOffspring("family1", 3)

retrieveGenotypes(24)
retrieveGenotypes(27)


