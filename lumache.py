"""
Lumache - Python library for cooks and food lovers.

This is a Python docstring, we can use reStructuredText syntax here!

.. code-block:: python

    # Import lumache
    import lumache

    # Call its only function
    get_random_ingredients(kind=["cheeses"])
"""

__version__ = "0.1.0"


class InvalidKindError(Exception):
    """Raised if the kind is invalid."""

    pass


def get_random_ingredients(kind=None):
    """
    Return a list of random ingredients as strings.

    :param kind: Optional "kind" of ingredients.
    :type kind: list[str] or None
    :raise lumache.InvalidKindError: If the kind is invalid.
    :return: The ingredients list.
    :rtype: list[str]
    """
    return ["shells", "gorgonzola", "parsley"]



class genome(object):
    '''
    Instantiates a genome object with genome-associated as paths pointing to the
    data source and set as class attributes. Data include the GFF3 file of 
    genome annotation (positional), InterproScan output of protein domains 
    identified on protein-coding genes (keyword), and directory with VCF files 
    of polymorphisms (keyword).
    
    :param gff3file: path to the GFF3 file
    :type gff3file: str
    :param iprfile: path to the InterproScan's output file
    :type iprfile: str
    :param vcffiles: path to the VCF files directory
    :type vcffiles: str
    '''

    _snpimpact = {'LOW':Color(0,150/255,50/255), 'MODERATE':Color(204/255,153/255,0/255),
                 'MODIFIER':Color(255/255,0/255,255/255),
                 'HIGH':Color(255/255,0,0), 'NOANN':Color(0,0,0)}


    def __init__(self, gff3file, iprfile=None, vcffiles=None):

        print('Initializing genome class...')
        genome.gff3file = gff3file
        genome.iprfile = iprfile
        genome.vcffiles = vcffiles
        logger.info(' genome class initialized')


    class gene(object):
        '''
        Instantiates a gene object with the method plot() to represent the 
        intron/exon structure of the gene from a GFF3 file, the protein domain 
        topology from InterproScan's output, and single nucleotide 
        polymofphisms (SNPs) from VCF files.

        :param mRNAid: gene identifier (ID) according to the GFF3 file annotations.
        :type mRNAid: str
        :param proteinid: protein identifier (ID) from the InterproScan output
        :type proteinid: str
        :param description: user-defined description of the gene
        :type description: str
        '''



        def __init__(self, mRNAid, proteinid=None, iprfile=None, description=None):

            print('Initializing gene class...')
            self.id = mRNAid
            self.proteinid = proteinid
            self.gff3file = genome.gff3file
            self.iprfile = genome.iprfile
            self.vcffiles = genome.vcffiles
            self.description = description
            self.db = gffutils.FeatureDB(genome.gff3file + '.db', keep_order=True); db = self.db

            assert db[mRNAid].featuretype == 'mRNA', 'only mRNA features are supported for \
                                                this version of the library'

            self.chrom = db[mRNAid].chrom
            self.start = db[mRNAid].start
            self.end = db[mRNAid].end
            logger.info(self.id + ' gene class initialized')


        def _transcriptpos_to_genomepos(self):
            '''
            Calculates genome coordinates for every nucleotide position
            of the transcript according to the GFF3 and FASTA
            files provided as input during the instantiation of the gene class.
            '''
            logger.info('starting...')
            db = self.db
            try:
                mrna = self.id
                coors = []
                lencdsaccum = 0
                for cds in db.children(db[mrna], featuretype='CDS', order_by='start'):
                    lencds = cds.end + 1 - cds.start
                    coors.append(list(range(cds.start, cds.end+1)))
                    lencdsaccum += lencds
                coors = sum(coors, [])
                if db[mrna].strand == '-':
                    coors.reverse()
                dcoors = {}
                for pos in range(lencdsaccum):
                    dcoors[pos+1] = coors[pos]
                self.coor = dcoors
                return dcoors
            except Exception as e:
                logger.error(e)
                raise

                
                
                
                
                
                
                
