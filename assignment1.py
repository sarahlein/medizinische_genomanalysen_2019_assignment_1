import mysql.connector
import os
import pybedtools
import pysam
import subprocess

__author__ = 'Sarah Meitz'

class Assignment1:
    
    def __init__(self):
        ## Your gene of interest
        self.target = "RRP1"
        self.gref = "hg38" #  gene reference,
        self.gene = self.download_gene_coordinates(self.gref, file_name="gene.txt") # hg38

        self.bamfile = os.path.join(os.getcwd(), "chr21.bam")
        self.baifile = os.path.join(os.getcwd(), "chr21.bam.bai")

        # Check if bam/bai file exist, if not download it
        # prior downloading of the bam file is highly recommended (large file!)
        if not os.path.isfile(self.bamfile):
            subprocess.call(["wget", "http://hmd.ait.ac.at/medgen2019/chr21.bam"])
        if not os.path.isfile(self.bamfile + ".bai"):
            subprocess.call(["samtools", "index", self.bamfile])

        self.samfile = pysam.AlignmentFile(self.bamfile, "rb")
        self.reads = list(self.samfile.fetch("chr21", self.gene.txStart, self.gene.txEnd))

    def download_gene_coordinates(self, genome_reference, file_name):

        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query for target gene
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)+ \
                " WHERE refGene.name2=" + '"' + self.target + '"' ""
        
        ## Execute query
        cursor.execute(query)

        ## Write to file
        with open(file_name, "w") as fh:
            for row in cursor:
                self.gene = gene(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8])
                fh.write(str(row) + "\n")

        # print statement for the gene class:
        #attribute_gene = vars(self.gene)
        #print(', '.join("%s: %s" % item for item in attribute_gene.items()))
            
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data")
        return self.gene
        
    def get_coordinates_of_gene(self):
        #get coordinates from UCSC file

        print("\nCoordinates of gene " + self.target + ": ")
        print("Start:\t", self.gene.txStart, "\nEnd:\t", self.gene.txEnd)
        print("todo")
        
    def get_gene_symbol(self):
        print("\nGene Symbol: ")
        print(self.gene.name)
        print()
                        
    def get_sam_header(self):
        print("Sam Header: ")
        for key, value in self.samfile.header['HD'].items():
            if key == "SO":
                print("Sorting order of alignments (SO): ", value)
            if key == "VN":
                print("Format version (VN): ", value)
            if key == "GO":
                print("Grouping of alignments (GO): ", value)
                print()
        
    def get_properly_paired_reads_of_gene(self):
        i = 0
        for read in self.reads:
            if read.is_proper_pair:
                i += 1
        if i == 0:
            print("\nNo properly paired end reads found \n")
        else:
            print("\nNumber of properly paired end reads:")
            print(i, "\n")
        
    def get_gene_reads_with_indels(self):
        i = 0
    # using cigar to find indels
        for read in self.reads:
            if not read.is_unmapped:
                cigar = read.cigar
                for (type, length) in cigar:

                    if (type == 1) or (type == 2):
                        i += 1
        if i == 0:
            print("No indels found! \n")
        else:
            print("Number of gene reads with indels:")
            print(i,"\n")
        
    def calculate_total_average_coverage(self):

        print("Start calculating the total average coverage, this may take a while...\n")

        a = pybedtools.BedTool(self.bamfile)
        b = a.genome_coverage(bg=True)

        i = 0
        average = 0

        for line in b:

            number = float(line[3])
            average += number
            i += 1

        coverage = average/i

        print("Total average coverage:")
        print("{0:.2f}".format(coverage), "\n")

    def calculate_gene_average_coverage(self):
        print("Start calculating the gene average coverage, this may take a while...\n")
        a = pybedtools.BedTool(self.bamfile)
        b = a.genome_coverage(bg=True)

        average = 0
        i = 0

        for line in b:
            number = float(line[3])
            cbeg = int(line[1])

            if cbeg > self.gene.txStart:
                if int(line[2]) <= self.gene.txEnd:
                    average += number
                    i += 1

        coverage = average / i

        print("Total gene average coverage:")
        print("{0:.2f}".format(coverage), "\n")
        
    def get_number_mapped_reads(self):

        i = 0
        for read in self.reads:

            if not read.is_unmapped:
                i += 1
        if i == 0:
            print("No mapped reads found \n")
        else:
            print("Number of mapped reads:")
            print(i, "\n")

    def get_region_of_gene(self):
        print("\nRegion of gene:")
        print("Chromosome: ", self.gene.chrom)
        
    def get_number_of_exons(self):

        if self.gene.exonCount == 0:
            print("\nNo exons found! \n")
        else:
            print("\nNumber of exons:")
            print(self.gene.exonCount, "\n")
    
    
    def print_summary(self):
        self.get_coordinates_of_gene()
        self.get_gene_symbol()
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_total_average_coverage()
        self.calculate_gene_average_coverage()
        self.get_number_mapped_reads()
        self.get_region_of_gene()
        self.get_number_of_exons()

        self.samfile.close()

# gene class stores all information about the target gene
class gene:
    def __init__(self, name2, name, chrom, txStart, txEnd, strand, exonCount, exonStarts, exonEnds):
        self.name2 = name2
        self.name = name
        self.chrom = chrom[3:]
        self.txStart = txStart
        self.txEnd = txEnd
        self.strand = strand
        self.exonCount = exonCount
        self.exonStarts = str(exonStarts).lstrip("b'").rstrip(",'").split(",")
        self.exonEnds = str(exonEnds).lstrip("b'").rstrip(",'").split(",")


def main():
    print("Assignment 1, by Sarah Meitz")
    assignment1 = Assignment1()
    assignment1.print_summary()
    
    
    print("Finished Assignment 1")
    
        
if __name__ == '__main__':
    main()
    
    
