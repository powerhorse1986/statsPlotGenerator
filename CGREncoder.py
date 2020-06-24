class GTFInfo:
    import pickle
    from os import path
    """
    parse the GTF file for retrieving information of each gene
    """
    def __init__(self, path_to_gtf):
        self.__revcomp = False # if the strand is positive, turn its to True
        self.__gtf_info = {} # keys are chromosomes, values are tuples (gene_id, int(start), int(stop + 1))

        if self.path.exists("gtf.pkl") == True and self.path.getsize("gtf.pkl") > 0:
            pkl_in = open("gtf.pkl", "rb")
            self.__gtf_info = self.pickle.load(pkl_in)
            pkl_in.close()
            print("GTF information loaded")
        else:
            gtf = open(path_to_gtf, "rb")
            for line in gtf:
                if "#" not in line:
                    chrom, src, seqtype, start, stop, score, strand, frame, attributes = line.strip().split("\t")

                    if "gene" in seqtype:
                        if strand == "-":
                            self.__revcomp = True
                        else:
                            self.__revcomp = False
                        if "gene_id" in attributes:
                            gene_id = attributes[9 : 24]
                            self.__gtf_info[gene_id] = (chrom, int(start) - 1, int(stop), self.__revcomp)

            print("GTF information loaded")
            pkl_out = open("gtf.pkl", 'wb')
            self.pickle.dump(self.__gtf_info, pkl_out)
            pkl_out.close()
            gtf.close()
            print("GTF information has been written to the local computer, it will be loaded directly next time")

    def getGtfInfo(self):
        return self.__gtf_info

class SeqOperator:
    from collections import defaultdict
    from os import path

    def __init__(self):
        self.__complementary_dict = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
        self.__key_set = set(self.__complementary_dict.keys())

    def reverse_seq(self, seq):
        seqList = list(seq)
        leftIdx = 0
        rightIdx = len(seqList) - 1

        while leftIdx <= rightIdx:
            tmpLeft = seqList[leftIdx]
            tmpRight = seqList[rightIdx]

            seqList[leftIdx] = tmpRight
            seqList[rightIdx] = tmpLeft

            leftIdx += 1
            rightIdx -= 1
        return "".join(seqList)

    def complement(self, seq):
        idx = 0
        seqList = list(seq)
        while idx < len(seqList):
            tmp = seqList[idx]
            tmp = self.__complementary_dict[tmp]
            seqList[idx] = tmp
            idx += 1
        return "".join(seqList)

    def complement_reverse(self, seq):
        leftIdx = 0
        rightIdx = len(seq) - 1
        seqList = list(seq)
        while leftIdx <= rightIdx:
            tmpLeft = seqList[leftIdx]
            tmpRight = seqList[rightIdx]
            if (tmpLeft not in self.__key_set) or (tmpRight not in self.__key_set):
                return "KeyError"
            tmpLeft = self.__complementary_dict[tmpLeft]
            tmpRight = self.__complementary_dict[tmpRight]

            seqList[leftIdx] = tmpRight
            seqList[rightIdx] = tmpLeft

            leftIdx += 1
            rightIdx -= 1
        return "".join(seqList)
    def countATGC(self, seq):
        """
        input:
        @seq: plain text sequence from loadSeq()

        return:
        a list of pertcentages of A, T, G, C and GC
        """
        dicA = defaultdict(int)
        dicT = defaultdict(int)
        dicG = defaultdict(int)
        dicC = defaultdict(int)
        seqLen = len(seq)
        for char in seq:
            if char == 'A':
                dicA['A'] = dicA['A'] + 1
            elif char == 'T':
                dicT['T'] = dicT['T'] + 1
            elif char == 'G':
                dicG['G'] = dicG['G'] + 1
            else:
                dicC['C'] = dicC['C'] + 1
        ratioA = dicA['A'] / seqLen
        ratioT = dicT['T'] / seqLen
        ratioG = dicG['G'] / seqLen
        ratioC = dicC['C'] / seqLen
        ratioGC = ratioG + ratioC
        print(ratioA)
        print(ratioT)
        print(ratioG)
        print(ratioC)
        print(ratioGC)
        return [ratioA, ratioT, ratioG, ratioC, ratioGC]

class SeqRetrieve:
    from pyensembl import EnsemblRelease
    from Bio import SeqIO
    from os import path
    import pickle
    import ensembl_rest as enr
    import sys
    import copy

    def __init__(self, ensembl_release = 75):
        self.__data = self.EnsemblRelease(ensembl_release)
        self.__transcript_info = None # __transcript_info contains transcripts informations, this is generated
                                      # by using pyensembl
        self.__gene = None
        self.__transcritSeqs = None # this __transcripts is a dictionary with keys are the ensembl transcripts ids
                                 # values are BioPython's SeqRecord objects
        self.__gene_coding_seq = None
        self.__gene_dict = {}
        self.__full_gene_seq = None
        self.__gtf_info = {}
        # self.__gtf_gene_info = None
        # self.__seqOperator = SeqOperator()

    def __remove_suf(self, record):
        """
        Given a sequence record, return the transcript id without the version suffix
        e.g.
        ENSG00000000419.1 -> ENSG00000000419
        """
        __full_id = record.id.split('.')
        #assert len(__full_id) == 2
        return __full_id[0]

    def getCodingSeq(self, geneID):
        if self.__transcritSeqs == None:
            self.__load_trans()
        self.__coding_seq_fetcher(geneID)
        return self.__gene_coding_seq

    def __load_trans(self):
        """
        __load_trans loads the .pkl file into a dictionary or parses the cdna.fa files into a dictionary
        """

        if self.path.exists('cdna.pkl') == True and self.path.getsize('cdna.pkl' > 0):
            pickle_in = open('cdna.pkl', 'rb')
            self.__transcritSeqs = self.pickle.load(pickle_in)
            pickle_in.close()
        else:
            self.__transcritSeqs = self.SeqIO.to_dict(SeqIO.parse('/home/lima/Project/simulation/cdna/Homo_sapiens.GRCh37.cdna.all.fa', 'fasta'),
                    key_function = self.__remove_suf)


    def __coding_seq_fetcher(self, geneID):
        """
        This fetcher loads a Genome object using ensemble geneID
        input:
        @geneID: a string, must be ensembl ID starts with ENSG, and no suffix
        return:
        a string contains all the transcripts sequences
        """
        self.__gene = self.__data.gene_by_id(geneID) # this is a Genome object contains only the basic information
                                                     # of this gene, includes exon, intro, start and end positions
                                                     # as well as the IDs of the transcripts contained in the gene

        self.__transcript_info = self.__gene.transcripts # Transcript object contains information of the transcripts
                                                         # contained in the gene
        __transcripts = self.__transcript_info

        self.__gene_coding_seq = ""
        for transcript in __transcripts:
            transcript_id = transcript.transcript_id
            tmp_seq = self.__transcritSeqs[transcript_id]
            self.__gene_coding_seq = self.__gene_coding_seq + tmp_seqs

    def getFullSeq(self, geneID):
        if len(self.__gene_dict) == 0:
            self.__load_gene_dict()
        if len(self.__gtf_info) == 0:
            self.__gtf_info = GTFInfo("/home/lima/Project/simulation/GRCh37/Homo_sapiens.GRCh37.87.gtf").getGtfInfo()
        self.__fetch_full_seq(geneID)
        return self.__full_gene_seq

    def __load_gene_dict(self):
        """
        loads the .pkl file into a dictionary or parses the GRCh37 master assembly fasta into a dictionary
        """
        if self.path.exists('dna_dict.pkl') == True and self.path.getsize('dna_dict.pkl' > 0):
            pkl_in = open('dna_dict.pkl', 'rb')
            self.__gene_dict = self.pickle.load(pkl_in)
            pkl_in.close()
        else:
            self.__gene_dict = self.SeqIO.to_dict(self.SeqIO.parse('/home/lima/Project/simulation/GRCh37/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa', 'fasta'))

    def __fetch_full_seq(self, geneID):
        __gtf_gene_info = self.__gtf_info[geneID]
        # print("gtf information length: ", len(__gtf_gene_info))
        # print(__gtf_gene_info)
        # assert len(__gtf_gene_info) == 4
        __gene_chrom = __gtf_gene_info[0]
        __gene_start = __gtf_gene_info[1]
        __gene_stop = __gtf_gene_info[2]
        __gene_rev = bool(__gtf_gene_info[3])
        tmp = str(self.__gene_dict[__gene_chrom].seq)
        tmp = tmp[__gene_start : __gene_stop]
        # print(type(tmp))
        # print(len(tmp))
        tmp = tmp.upper()
        # cr_processed_tmp = ""
        cr_processed_tmp = SeqOperator().complement_reverse(tmp)
        self.__full_gene_seq = cr_processed_tmp

    def debug(self, geneID):
        if len(self.__gene_dict) == 0:
            self.__load_gene_dict()
        if len(self.__gtf_info) == 0:
            self.__gtf_info = GTFInfo("/home/lima/Project/simulation/GRCh37/Homo_sapiens.GRCh37.87.gtf").getGtfInfo()
        self.__fetch_debug_seq(geneID)
        return self.__full_gene_seq

    def __fetch_debug_seq(self, geneID):
        __gtf_gene_info = self.__gtf_info[geneID]
        # print("gtf information length: ", len(__gtf_gene_info))
        # print(__gtf_gene_info)
        # assert len(__gtf_gene_info) == 4
        __gene_chrom = __gtf_gene_info[0]
        __gene_start = __gtf_gene_info[1]
        __gene_stop = __gtf_gene_info[2]
        # __gene_rev = bool(__gtf_gene_info[3])
        tmp = str(self.__gene_dict[__gene_chrom].seq)
        tmp = tmp[__gene_start : __gene_stop]
        # print(type(tmp))
        # print(len(tmp))
        tmp = tmp.upper()
        self.__full_gene_seq = tmp




class SeqLoader:
    import requests

    """
    # input:
    # @GRCh37: bool
    """
    def __init__(self, GRCh37 = True):
        if GRCh37 == True:
            self.__server = "http://grch37.rest.ensembl.org"
            self.__ext = ""
            self.__seq = ""
            # self.__geneID = geneID

    def __loadSeq(self, geneID):
        self.__ext = "/sequence/id/{gene_id}?".format(gene_id = geneID)
        seq = self.requests.get(self.__server + self.__ext, headers={ 'content-Type' : 'text/plain'})
        if not seq.ok:
          seq.raise_for_status()
          sys.exit()
        self.__seq = seq.text

    def getSeq(self, geneID):
        self.__loadSeq(geneID)
        return self.__seq

class CGREncoder:
    import pandas as pd
    import math
    import numpy as np
    from collections import defaultdict
    def __init__(self):
        self.__kmer_counts = None
        self.__kmer_prob = None
        self.__kmerSize = -1
        self.__matrixsize = -1
        self.__chaos = None

    def __kmerGenerator(self, seq, k):
        """
        input:
        @seq: string, a plain text sequence
        @k: a positive integer

        return:
        self.__kmer: a dictionary with keys being the kmers, values being the counts of kmers
        """
        self.__kmer_counts = self.defaultdict(int)
        self.__kmer_prob = self.defaultdict(int)
        # print("Generating kmer pool")
        size = int(self.math.sqrt(4 ** k))
        self.__matrixsize = size
        self.__chaos = self.np.zeros((size, size), dtype = self.np.float)
        self.__kmerSize = len(seq) - k + 1
        for i in range(self.__kmerSize):
            self.__kmer_counts[seq[i : i + k]] += 1

        # some sequences do have Ns, eventhough only one
        for key in list(self.__kmer_counts.keys()):
            if 'N' in key:
                del self.__kmer_counts[key]

        for key, value in self.__kmer_counts.items():
            self.__kmer_prob[key] = float(value) / (len(seq) - k + 1)

    def getKmerSize(self):
        return self.__kmerSize

    def getKmer(self):
        return [self.__kmer_counts, self.__kmer_prob]

    def encoding(self, seq, k):
        self.__kmerGenerator(seq, k)
        for key, value in self.__kmer_prob.items():
            minx = 0
            miny = 0
            maxx = self.__matrixsize - 1
            maxy = self.__matrixsize - 1
            midx = self.__midPoint(minx, maxx)
            midy = self.__midPoint(miny, maxy)
            charIdx = 0
            self.__helper(charIdx, minx, midx, maxx, miny, midy, maxy, key, value)

    def getChaosMatrix(self):
        return self.__chaos.copy()

    def __midPoint(self, small, large):
        return int(small + (large - small) / 2)

    def __helper(self, charIdx, minx, midx, maxx, miny, midy, maxy, key, value):
        if ((minx == maxx or miny == maxy) or charIdx >= len(key)):
            # idxx = self.__midPoint(minx, maxx)
            # idxy = self.__midPoint(miny, maxy)
            # print(idxy)
            # print(idxx)
            self.__chaos[minx][miny] = value
            return
        #print(key)
        char = key[charIdx]
        #print(char)
        charIdx = charIdx + 1
        if char == 'A':
            minx = minx
            maxx += midx
            miny = miny
            maxy += midy

        elif char == 'T':
            minx += midx
            maxx = maxx
            miny = miny
            maxy += midy

        elif char == 'G':
            minx += midx
            maxx = maxx
            miny += midy
            maxy = maxy
        else:
            minx = minx
            maxx += midx
            miny += midy
            maxy = maxy
        midx = self.__midPoint(minx, maxx)
        midy = self.__midPoint(miny, maxy)

        self.__helper(charIdx, minx, midx, maxx, miny, midy, maxy, key, value)

class EncodingByTreatMent:
    import pandas as pd
    import multiprocessing as mp
    import KVRPlot
    __trtmnt = {'high': pd.read_csv('/home/lima/Project/simulation/Comparision/high.txt', skiprows = 0).values.tolist(),
                     'low': pd.read_csv('/home/lima/Project/simulation/Comparision/low.txt', skiprows = 0).values.tolist(),
                     'veh': pd.read_csv('/home/lima/Project/simulation/Comparision/veh.txt', skiprows = 0).values.tolist(),
                     'pt': pd.read_csv('/home/lima/Project/simulation/Comparision/pt.txt', skiprows = 0).values.tolist(),
                     'all': pd.read_csv('/home/lima/Project/simulation/Comparision/samples.txt', skiprows = 0).values.tolist()}

    def __init__(self, trtmnt, normalized, pre_process_software, cut):
        self.__loader = self.KVRPlot.DataLoader(trtmnt, normalized)
        # self.__kall = self.__loader.getKall()
        if pre_process_software == "RSEM":
            __temp = self.__loader.getRSEM()
            self.__counts_matrix = __temp[cut]
            # print(self.__rsem)
        else:
            self.__counts_matrix = self.__loader.getKall()
            
        self.__genes = list(self.__counts_matrix.index.values)
        self.__gene_number = len(self.__genes)
        # print(len(self.__genes))
        # self.__cut = cut
        self.__samples = self.__trtmnt[trtmnt]
        # print(self.__samples)
        self.__encoder = CGREncoder()
        # self.__loader = SeqLoader()
        self.__retriever = SeqRetrieve()
        self.__trtmnt_grp = trtmnt
        self.__result = self.mp.Manager().dict()

    def convention_convert(self, k):
        for sample in self.__samples:
            __tmpDF = self.__counts_matrix[sample]
            __tmpDict = {}
            for index in range(len(__tmpDF)):
                __gene = self.__genes[index]
                print("Working with {gene}".format(gene = __gene))
                #__tmpSeq = self.__loader.getSeq(__gene)
                self.__encoder.encoding(__gene, k)
                __tmpMatrix = self.__encoder.getChaosMatrix()
                __tmpMatrix = __tmpMatrix.reshape(1, 4 ** k)
                __tmpList = __tmpMatrix.tolist()
                __tmpDict[__gene]: __tmpList
            self.__result[sample]: __tmpDict

    def __multicoreConvert(self, k, index):
        __gene = self.__genes[index]
        print("Working with {gene}".format(gene = __gene))
        __tmpSeq = self.__retriever.getFullSeq(__gene)
        self.__encoder.encoding(__tmpSeq, k)
        __tmpMatrix = self.__encoder.getChaosMatrix()
        __tmpMatrix = __tmpMatrix.reshape(1, 4 ** k)
        __tmpList = tuple(__tmpMatrix.tolist())
        # __tmpDict[__gene] = __tmpList
        self.__result[__gene] = __tmpList
        # result = {__gene: __tmpList}
        # return result

    def multicore_convert(self, k, multiprocessing):
        pool = self.mp.Pool(processes = multiprocessing)
        # for sample in self.__samples:
        # __tmpDF = self.__rsem[sample]
        for x in range(self.__gene_number):
            pool.apply_async(self.__multicoreConvert(2, x), args = (x, ))
        # results = tuple(results)
        # output = []
        # for result in results:
        #    output.append(result.get())
        pool.close()
        pool.join()
        # output = [p.get() for p in results]
        # self.__result[self.__trtmnt_grp] = self.__result
        print("{Group} being converted.".format(Group = self.__trtmnt_grp))

    def getConvertedGroup(self):
        return self.__result
