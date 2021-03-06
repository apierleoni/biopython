# Copyright 2010 by Andrea Pierleoni
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "uniprot" XML file formats.

See also:

http://www.uniprot.org
"""
import sys

from Bio import Seq
from Bio import SeqFeature
from Bio import Alphabet
from Bio.SeqRecord import SeqRecord
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings
try:
    if (3,0,0) <= sys.version_info < (3,1,2):
        #workaround for bug in python 3 to 3.1.2  see http://bugs.python.org/issue9257
        from xml.etree import ElementTree as ElementTree
    else:
        from xml.etree import cElementTree as ElementTree
except ImportError:
    try:
        from xml.etree import ElementTree as ElementTree
    except ImportError:
        # Python 2.4 -- check for 3rd-party implementations
        try:
            from lxml import etree as ElementTree
        except ImportError:
            try:
                import cElementTree as ElementTree
            except ImportError:
                try:
                    from elementtree import ElementTree
                except ImportError:
                    from Bio import MissingExternalDependencyError
                    raise MissingExternalDependencyError(
                            "No ElementTree module was found. "
                            "Use Python 2.5+, lxml or elementtree if you "
                            "want to use Bio.SeqIO.UniprotIO.")

NS = "{http://uniprot.org/uniprot}"
REFERENCE_JOURNAL = "%(name)s %(volume)s:%(first)s-%(last)s(%(pub_date)s)"

def UniprotIterator(handle, alphabet=Alphabet.ProteinAlphabet(), return_raw_comments=False, skip_parsing_errors=False):
    '''Generator Function
    parses an XML entry at a time from any UniProt XML file 
    returns a SeqRecord for each iteration
    
    This generator can be used in Bio.SeqIO
    
    return_raw_comments = True --> comment fields are returned as complete xml to allow further processing
    skip_parsing_errors = True --> if parsing errors are found, skip to next entry
    '''
    if isinstance(alphabet, Alphabet.NucleotideAlphabet):
        raise ValueError, "Wrong alphabet %r" % alphabet
    if isinstance(alphabet, Alphabet.Gapped):
        if isinstance(alphabet.alphabet, Alphabet.NucleotideAlphabet):
            raise ValueError, "Wrong alphabet %r" % alphabet

    if not hasattr(handle, "read"):
        if type(handle)==type(''):
            handle=StringIO(handle)
        else:
            raise Exception('An XML-containing handler or an XML string must be passed')

    for event, elem in ElementTree.iterparse(handle, events=("start", "end")):
        if event=="end" and elem.tag == NS + "entry":
            try:
                yield Parser(elem, alphabet=alphabet, return_raw_comments=return_raw_comments).parse()
                elem.clear()
            except:
                # TODO: log warning, test
                if skip_parsing_errors:
                    pass
                else:
                    raise

class Parser():
    '''Parse a UniProt XML entry to a SeqRecord
    return_raw_comments=True to get back the complete comment field in XML format
    alphabet=Alphabet.ProteinAlphabet()    can be modified if needed, default is protein alphabet.
    '''
    def __init__(self, elem, alphabet=Alphabet.ProteinAlphabet(), return_raw_comments=False):
        self.entry=elem
        self.alphabet=alphabet
        self.return_raw_comments=return_raw_comments
    
    def parse(self):
        '''parse the input '''
        assert self.entry.tag == NS + 'entry'
        
        def append_to_annotations(key, value):
            if not self.ParsedSeqRecord.annotations.has_key(key):
                self.ParsedSeqRecord.annotations[key]=[]
            if value not in self.ParsedSeqRecord.annotations[key]:
                self.ParsedSeqRecord.annotations[key].append(value)
            
        def _parse_name(element):
            '''use name as name'''
            self.ParsedSeqRecord.name=element.text
            '''add name to dbxrefs'''
            self.ParsedSeqRecord.dbxrefs.append(self.dbname+':'+element.text)
        
        def _parse_accession(element):
            append_to_annotations('accessions', element.text)# to cope with SwissProt plain text parser
            '''add accessions to dbxrefs'''
            self.ParsedSeqRecord.dbxrefs.append(self.dbname+':'+element.text)
        
        def _parse_protein(element):
            '''Parse protein names'''
            descr_set=False
            for protein_element in element.getchildren():
                if protein_element.tag in [NS + 'recommendedName', NS + 'alternativeName']:#recommendedName tag are parsed before
                    '''use protein fields for name and description '''
                    for rec_name in protein_element.getchildren():
                        ann_key='%s_%s' % (protein_element.tag.replace(NS,''), rec_name.tag.replace(NS,''))
                        append_to_annotations(ann_key, rec_name.text)
                        if (rec_name.tag==NS + 'fullName') and not descr_set:
                            self.ParsedSeqRecord.description=rec_name.text
                            descr_set=True
                elif protein_element.tag==NS + 'component':
                    pass #not parsed 
                elif protein_element.tag==NS + 'domain':
                    pass #not parsed 
        
        def _parse_gene(element):
            for genename_element in element.getchildren():  
                if genename_element.attrib.has_key('type'):
                    ann_key='gene_%s_%s' % (genename_element.tag.replace(NS,''), genename_element.attrib['type'])
                    if genename_element.attrib['type']=='primary':
                        self.ParsedSeqRecord.annotations[ann_key]=genename_element.text
                    else:
                        append_to_annotations(ann_key,genename_element.text)
        
        def _parse_geneLocation(element):
            append_to_annotations('geneLocation', element.attrib['type'])
        
        def _parse_organism(element):
            com_name=sci_name=''
            for organism_element in element.getchildren():  
                if organism_element.tag==NS + 'name':
                    if organism_element.attrib['type']== 'scientific':
                        sci_name=organism_element.text
                    elif organism_element.attrib['type']== 'common':
                        com_name=organism_element.text
                    else:
                        append_to_annotations("organism_name", organism_element.text)
                elif organism_element.tag==NS + 'dbReference':
                    self.ParsedSeqRecord.dbxrefs.append(organism_element.attrib['type']+':'+organism_element.attrib['id'])
                elif organism_element.tag==NS + 'lineage':
                    for taxon_element in organism_element.getchildren():
                        if taxon_element.tag==NS + 'taxon':
                            append_to_annotations('taxonomy',taxon_element.text)
            organism_name=sci_name 
            if com_name:
                organism_name+=' (%s)'%com_name
            self.ParsedSeqRecord.annotations['organism']=organism_name
            
        def _parse_organismHost(element):
            for organism_element in element.getchildren():  
                if organism_element.tag==NS + 'name': 
                    append_to_annotations("organismHost_name", organism_element.text)
                        
        def _parse_keyword(element):      
            append_to_annotations('keywords',element.text)
        
        def _parse_comment(element):
            '''Comment fields are very heterogeneus. each type has his own (frequently mutated) schema.
            To store all the contained data, more complex data structures are needed, such as 
            annidated dictionaries. This is left to end user, by optionally setting:
            
            return_raw_comments=True 
            
            the orginal XMLs is returned in the annotation fields.
            
            available comment types at december 2009:
                "allergen"
                "alternative products"
                "biotechnology"
                "biophysicochemical properties"
                "catalytic activity"
                "caution"
                "cofactor"
                "developmental stage"
                "disease"
                "domain"
                "disruption phenotype"
                "enzyme regulation"
                "function"
                "induction"
                "miscellaneous"
                "pathway"
                "pharmaceutical"
                "polymorphism"
                "PTM"
                "RNA editing"
                "similarity"
                "subcellular location"
                "sequence caution"
                "subunit"
                "tissue specificity"
                "toxic dose"
                "online information"
                "mass spectrometry"
                "interaction"
            '''
            
            simple_comments=["allergen",
                            "biotechnology",
                            "biophysicochemical properties",
                            "catalytic activity",
                            "caution",
                            "cofactor",
                            "developmental stage",
                            "disease",
                            "domain",
                            "disruption phenotype",
                            "enzyme regulation",
                            "function",
                            "induction",
                            "miscellaneous",
                            "pathway",
                            "pharmaceutical",
                            "polymorphism",
                            "PTM",
                            "RNA editing",#positions not parsed
                            "similarity",
                            "subunit",
                            "tissue specificity",
                            "toxic dose",
                             ]

            if element.attrib['type'] in simple_comments:
                ann_key='comment_%s' % element.attrib['type'].replace(' ','')
                for text_element in element.getiterator(NS + 'text'):
                    if text_element.text:
                        append_to_annotations(ann_key,text_element.text)
            elif element.attrib['type']=='subcellular location':
                for subloc_element in element.getiterator(NS + 'subcellularLocation'):
                    for el in subloc_element.getchildren():
                        if el.text:
                            ann_key='comment_%s_%s' % (element.attrib['type'].replace(' ',''), el.tag.replace(NS,''))
                            append_to_annotations(ann_key,el.text)
            elif element.attrib['type']=='interaction':
                for interact_element in element.getiterator(NS +'interactant'):
                    ann_key='comment_%s_intactId' % element.attrib['type']
                    append_to_annotations(ann_key,interact_element.attrib['intactId'])
            elif element.attrib['type']=='alternative products':
                for alt_element in element.getiterator(NS +'isoform'):
                    ann_key='comment_%s_isoform' % element.attrib['type'].replace(' ','')
                    for id_element in alt_element.getiterator(NS +'id'):
                        append_to_annotations(ann_key,id_element.text)
            elif element.attrib['type']=='mass spectrometry':
                ann_key='comment_%s' % element.attrib['type'].replace(' ','')
                start=end=0
                for loc_element in element.getiterator(NS +'location'):
                    pos_els=loc_element.getiterator(NS +'position')
                    pos_els=list(pos_els)
                    # this try should be avoided, maybe it is safer to skip postion parsing for mass spectrometry
                    try:
                        if pos_els:
                            start=end=int(pos_els[0].attrib['position'])-1
                        else:
                                start=int(loc_element.getiterator(NS +'begin')[0].attrib['position'])-1
                                end=int(loc_element.getiterator(NS +'end')[0].attrib['position'])-1
                    except :#undefined positions or erroneusly mapped
                        pass    
                mass=element.attrib['mass']
                method=element.attrib['mass']
                if start==end==0:  
                    append_to_annotations(ann_key,'undefined:%s|%s'%(mass,method))
                else:
                    append_to_annotations(ann_key,'%s..%s:%s|%s'%(start,end,mass,method))
            elif element.attrib['type']=='sequence caution':
                pass#not parsed: few information, complex structure
            elif element.attrib['type']=='online information':
                for link_element in element.getiterator(NS +'link'):
                    ann_key='comment_%s' % element.attrib['type'].replace(' ','')
                    for id_element in link_element.getiterator(NS +'link'):
                        append_to_annotations(ann_key,'%s@%s'%(element.attrib['name'],link_element.attrib['uri']))            
            
            '''return raw XML comments if needed '''
            if self.return_raw_comments:
                ann_key='comment_%s_xml' % element.attrib['type'].replace(' ','')
                append_to_annotations(ann_key,ElementTree.tostring(element))
                
        
        def _parse_dbReference(element):
            self.ParsedSeqRecord.dbxrefs.append(element.attrib['type']+':'+element.attrib['id'])
            '''<dbReference type="PDB" key="11" id="2GEZ">
               <property value="X-ray" type="method"/>
               <property value="2.60 A" type="resolution"/>
               <property value="A/C/E/G=1-192, B/D/F/H=193-325" type="chains"/>
             </dbReference>'''
            if 'type' in element.attrib:
                if element.attrib['type'] == 'PDB':
                        method=""
                        resolution=""
                        for ref_element in element.getchildren():  
                            if ref_element.tag==NS + 'property':
                                dat_type=ref_element.attrib['type']
                                if dat_type=='method':
                                    method=ref_element.attrib['value']
                                if dat_type=='resolution':
                                    resolution=ref_element.attrib['value']
                                if dat_type=='chains':
                                    pairs=ref_element.attrib['value'].split(',')
                                    for elem in pairs:
                                        pair=elem.strip().split('=')
                                        if pair[1]!='-':
                                            feature=SeqFeature.SeqFeature()
                                            feature.type=element.attrib['type']
                                            feature.qualifiers['name']=element.attrib['id']
                                            feature.qualifiers['method']=method
                                            feature.qualifiers['resolution']=resolution
                                            feature.qualifiers['chains']=pair[0].split('/')
                                            start=int(pair[1].split('-')[0])-1
                                            end=int(pair[1].split('-')[1])-1
                                            feature.location=SeqFeature.FeatureLocation(start,end)
                                            self.ParsedSeqRecord.features.append(feature)

            for ref_element in  element.getchildren():  
                if ref_element.tag==NS + 'property':
                    pass# this data cannot be fitted in a seqrecord object with a simple list. however at least ensembl and EMBL parsing can be improved to add entries in dbxrefs
            
        def _parse_reference(element):
            reference=SeqFeature.Reference()
            authors=[]
            scopes=[]
            tissues=[]
            journal_name=''
            pub_type=''
            pub_date=''
            for ref_element in element.getchildren():
                if ref_element.tag==NS + 'citation':
                    pub_type=ref_element.attrib['type']
                    if pub_type=='submission':
                        pub_type+=' to the '+ref_element.attrib['db']
                    if ref_element.attrib.has_key('name'):
                        journal_name=ref_element.attrib['name']
                    if ref_element.attrib.has_key('date'):
                        pub_date=ref_element.attrib['date']
                    else:
                        pub_date=''
                    if ref_element.attrib.has_key('volume'):
                        j_volume=ref_element.attrib['volume']
                    else:
                        j_volume=''
                    if ref_element.attrib.has_key('first'):
                        j_first=ref_element.attrib['first']
                    else:
                        j_first=''
                    if ref_element.attrib.has_key('last'):
                        j_last=ref_element.attrib['last']
                    else:
                        j_last=''
                    for cit_element in ref_element.getchildren():
                        if cit_element.tag==NS + 'title':
                            reference.title=cit_element.text
                        elif cit_element.tag==NS + 'authorList':
                            for person_element in cit_element.getchildren():
                                authors.append(person_element.attrib['name'])
                        elif cit_element.tag==NS + 'dbReference':
                            self.ParsedSeqRecord.dbxrefs.append(cit_element.attrib['type']+':'+cit_element.attrib['id'])
                            if cit_element.attrib['type']=='PubMed':
                                reference.pubmed_id=cit_element.attrib['id']
                            elif ref_element.attrib['type']=='MEDLINE':
                                reference.medline_id=cit_element.attrib['id']
                elif ref_element.tag==NS + 'scope':
                    scopes.append(ref_element.text)
                elif ref_element.tag==NS + 'source':
                    for source_element in ref_element.getchildren():
                        if source_element.tag==NS + 'tissue':
                            tissues.append(source_element.text)
            if scopes:
                scopes_str='Scope: '+', '.join(scopes)
            else:
                scopes_str=''
            if tissues:
                tissues_str='Tissue: '+', '.join(tissues)
            else:
                tissues_str=''
            
            reference.location = [] #locations cannot be parsed since they are actually written in free text inside scopes so all the references are put in the annotation.
            reference.authors = ', '.join(authors) 
            if journal_name:
                if pub_date and j_volume and j_first and j_last:
                    reference.journal = REFERENCE_JOURNAL % dict(name=journal_name,
                        volume=j_volume, first=j_first, last=j_last, pub_date=pub_date)
                else:
                    reference.journal = journal_name 
            reference.comment = ' | '.join((pub_type,pub_date,scopes_str,tissues_str))
            append_to_annotations('references', reference)
            
        def _parse_feature(element):
            feature=SeqFeature.SeqFeature()
            for k,v in element.attrib.items():
                feature.qualifiers[k]=v
            if element.attrib.has_key('type'):
                feature.type=element.attrib['type']
            else:
                feature.type=''
            if element.attrib.has_key('type'):
                feature.type=element.attrib['type']
            if element.attrib.has_key('id'):
                feature.id=element.attrib['id']
            for feature_element in element.getchildren():
                if feature_element.tag==NS + 'location':
                    try:
                        position_elements=feature_element.findall(NS + 'position')
                        if position_elements:
                            position=int(position_elements[0].attrib['position'])-1
                            feature.location=SeqFeature.FeatureLocation(position,position)
                        else:
                            start_positions_elements=feature_element.findall(NS + 'begin')
                            
                            if start_positions_elements:
                                if start_positions_elements[0].attrib.has_key('position'):
                                    start_position=int(start_positions_elements[0].attrib['position'])-1
                                elif start_positions_elements[0].attrib.has_key('status'):#fuzzy location
                                    return #skip feature with unparsable position
                            end_positions_elements=feature_element.findall(NS + 'end')
                            if end_positions_elements:
                                if end_positions_elements[0].attrib.has_key('position'):
                                    end_position=int(end_positions_elements[0].attrib['position'])-1
                                elif end_positions_elements[0].attrib.has_key('status'):#fuzzy location
                                    return #skip feature with unparsable position
                            feature.location=SeqFeature.FeatureLocation(start_position,end_position)
                    except:
                        return #skip feature with unparsable position
                else:
                    try:
                        feature.qualifiers[feature_element.tag.replace(NS,'')]=feature_element.text
                    except:
                        pass#skip unparsable tag
            self.ParsedSeqRecord.features.append(feature)
            
        def _parse_proteinExistence(element):
            append_to_annotations('proteinExistence', element.attrib['type'])   
            
        def _parse_evidence(element):
            for k, v in  element.attrib.items():
                ann_key = k
                append_to_annotations(ann_key, v)   
        
        def  _parse_sequence(element):
            for k, v in element.attrib.items():
                if k in ("length", "mass", "version"):
                    self.ParsedSeqRecord.annotations['sequence_%s' % k] = int(v)
                else:
                    self.ParsedSeqRecord.annotations['sequence_%s' % k] = v
            seq=''.join((element.text.split()))
            self.ParsedSeqRecord.seq=Seq.Seq(seq,self.alphabet)
            
        #============================================#
        '''Initialize SeqRecord '''
        self.ParsedSeqRecord=SeqRecord('', id='') 
        
        '''Entry attribs parsing '''
        if self.entry.attrib.has_key('dataset'):
            self.dbname=self.entry.attrib['dataset']
        else:
            self.dbname='UnknownDataset'#this should not happen!
        '''add attribs to annotations '''
        for k, v in self.entry.attrib.items():
            if k in ("version"):
                '''original'''
                #self.ParsedSeqRecord.annotations["entry_%s" % k] = int(v)
                '''to cope with swissProt plain text parser. this can cause errors 
                if the attrib has the same name of an other annotation'''
                self.ParsedSeqRecord.annotations[k] = int(v)
            else:
                #self.ParsedSeqRecord.annotations["entry_%s" % k] = v
                self.ParsedSeqRecord.annotations[k] = v # to cope with swissProt plain text parser

        '''Top-to-bottom entry children parsing '''
        for element in self.entry.getchildren():
            if element.tag==NS + 'name':
                _parse_name(element)
            elif element.tag==NS + 'accession':
                _parse_accession(element)
            elif element.tag==NS + 'protein':
                _parse_protein(element)  
            elif element.tag==NS + 'gene':
                _parse_gene(element)
            elif element.tag==NS + 'geneLocation':
                _parse_geneLocation(element)
            elif element.tag==NS + 'organism':
                _parse_organism(element)          
            elif element.tag==NS + 'organismHost':
                _parse_organismHost(element)
            elif element.tag==NS + 'keyword':
                _parse_keyword(element)
            elif element.tag==NS + 'comment':
                _parse_comment(element)
            elif element.tag==NS + 'dbReference':
                _parse_dbReference(element)
            elif element.tag==NS + 'reference':
                _parse_reference(element)
            elif element.tag==NS + 'feature':
                _parse_feature(element)
            elif element.tag==NS + 'proteinExistence':
                _parse_proteinExistence(element)
            elif element.tag==NS + 'evidence':
                _parse_evidence(element)
            elif element.tag==NS + 'sequence':
                _parse_sequence(element)
            else:
                pass   
            
        self.ParsedSeqRecord.dbxrefs=list(set(self.ParsedSeqRecord.dbxrefs))#remove duplicate dbxrefs
        self.ParsedSeqRecord.dbxrefs.sort()

        # use first accession as id
        if not self.ParsedSeqRecord.id:
            self.ParsedSeqRecord.id=self.ParsedSeqRecord.annotations['accessions'][0]


        
        return self.ParsedSeqRecord
        
