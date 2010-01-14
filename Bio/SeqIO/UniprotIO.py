'''
Created on 04/dic/2009

@author: andreapierleoni


Bio.UniprotIO support for the "uniprot" (aka UniProt XML) file format.

You are expected to use this module via the Bio.SeqIO functions.


'''

from Bio import Seq
from Bio import SeqFeature
from Bio import Alphabet
from Bio.SeqRecord import SeqRecord
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
from xml.etree import ElementTree
import warnings


def UniprotIterator(handle,root_element='entry',alphabet=Alphabet.ProteinAlphabet(),return_raw_comments=False, skip_parsing_errors=False):
    '''Generator Function
    parses an XML entry at a time from any UniProt XML file 
    returns a SeqRecord for each iteration
    
    This generator can be used in Bio.SeqIO
    
    return_raw_comments = True --> comment fields are returned as complete xml to allow further processing
    skip_parsing_errors = True --> if parsing errors are found, skip to next entry
    
    '''
    if not hasattr(handle, "read"):
        if type(handle)==type(''):
            handle=StringIO(handle)
        else:
            raise Exception('An XML-containing handler or an XML string must be passed')
    single_entry=[]
    write=False#this check is needed to save memory by loading only data encolsed between root_element tags. if a malformed file or a wrong format file is passed, it will not be loaded completely in memory
    for line in handle:
        if '<%s'%root_element in line:
            single_entry=[]
            write=True
        if write:
            single_entry.append(line)
        if '</%s>'%root_element in line:
            write=False
            try:
                yield Parser(''.join(single_entry),alphabet=alphabet,return_raw_comments=return_raw_comments).parse()
            except Exception,error:
                if skip_parsing_errors:
                    warnings.warn('Error in parsing xml format: %s'%error)
                else:
                    raise ValueError('Error in parsing xml format: %s'%error)

    return 


class Parser():
    '''parse a UniProt XML entry to a SeqRecord
    return_raw_comments=True to get back the complete comment field in XML format
    alphabet=ProteinAlphabet()    can be modified if needed, default is protein alphabet.
    '''
    
    def __init__(self,xml_entry,alphabet=Alphabet.ProteinAlphabet(),return_raw_comments=False):
        
        if not hasattr(xml_entry, "read"):
            xml_entry=StringIO(xml_entry)
        self.raw_xml=xml_entry
        self.alphabet=alphabet
        self.return_raw_comments=return_raw_comments
        
        self.validate()
    
    def validate(self):
        try:
            self.dom=ElementTree.parse(self.raw_xml)
            self.entry=self.dom.getroot()
            if self.entry.tag!='entry':
                raise Exception('Malformed XML, the root should be an "entry" tag')
            self.validated=True
        except:
            self.validated=False

        
    def parse(self):
        
        if not self.validated:
            return None
        
        
        self.ParsedSeqRecord=SeqRecord('')    
        
        def append_to_annotations(key,value):
            if not self.ParsedSeqRecord.annotations.has_key(key):
                self.ParsedSeqRecord.annotations[key]=[]
            if value not in self.ParsedSeqRecord.annotations[key]:
                self.ParsedSeqRecord.annotations[key].append(value)
            
        def _parse_name(element):
            '''use name as id'''
            self.ParsedSeqRecord.id=element.text
            '''add name to dbxrefs'''
            self.ParsedSeqRecord.dbxrefs.append(self.dbname+':'+element.text)
        
        def _parse_accession(element):
            '''add accessions to dbxrefs'''
            self.ParsedSeqRecord.dbxrefs.append(self.dbname+':'+element.text)
        
        def _parse_protein(element):
            '''Parse protein names'''
            descr_set=False
            name_set=False
            for protein_element in element.getchildren():
                if protein_element.tag in ['recommendedName','alternativeName']:#recommendedName tag are parsed before
                    '''use protein fields for name and description '''
                    for rec_name in  protein_element.getchildren():
                        ann_key='_'.join((protein_element.tag,rec_name.tag))
                        append_to_annotations(ann_key,rec_name.text)
                        if (rec_name.tag=='fullName') and not descr_set:
                            self.ParsedSeqRecord.description=rec_name.text
                            descr_set=True
                        elif (rec_name.tag=='shortName') and not name_set:
                            self.ParsedSeqRecord.name=rec_name.text
                            name_set=True
                elif protein_element.tag=='component':
                    pass #not parsed 
                elif protein_element.tag=='domain':
                    pass #not parsed 
        
        def _parse_gene(element):
            for genename_element in element.getchildren():  
                if genename_element.attrib.has_key('type'):
                    ann_key='_'.join((element.tag,genename_element.tag,genename_element.attrib['type']))
                    if genename_element.attrib['type']=='primary':
                        self.ParsedSeqRecord.annotations[ann_key]=genename_element.text
                    else:
                        append_to_annotations(ann_key,genename_element.text)
        
        def _parse_geneLocation(element):
            append_to_annotations(element.tag,element.attrib['type'])
        
        def _parse_organism(element):
            com_name=sci_name=''
            for organism_element in element.getchildren():  
                if organism_element.tag=='name':
                    if organism_element.attrib['type']== 'scientific':
                        sci_name=organism_element.text
                    elif organism_element.attrib['type']== 'common':
                        com_name=organism_element.text
                    else:
                        ann_key='_'.join((element.tag,organism_element.tag))
                        append_to_annotations(ann_key,organism_element.text)
                elif organism_element.tag=='dbReference':
                    self.ParsedSeqRecord.dbxrefs.append(organism_element.attrib['type']+':'+organism_element.attrib['id'])
                elif organism_element.tag=='lineage':
                    for taxon_element in organism_element.getchildren():
                        if taxon_element.tag=='taxon':
                            append_to_annotations('taxonomy',taxon_element.text)
            organism_name=sci_name 
            if com_name:
                organism_name+=' (%s)'%com_name
            self.ParsedSeqRecord.annotations['organism']=organism_name
            
        def _parse_organismHost(element):
            for organism_element in element.getchildren():  
                if organism_element.tag=='name': 
                    ann_key='_'.join((element.tag,organism_element.tag))
                    append_to_annotations(ann_key,organism_element.text)
                        
        def _parse_keyword(element):      
            ann_key='keywords'
            append_to_annotations(ann_key,element.text)
        
        def _parse_comment(element):
            '''Comment fields are very heterogeneus. each type has his own (frequently mutated) schema.
            To store all the contained data, more complex data structures are needed, such as 
            annidated dictionaires, but this is ledft to to end user.
            By optionally setting:
            
            return_raw_comments=True 
            
            the orginal XMLs is returned in the annotation fields.
            Note that this will behave badly when saving the seqrecord as genbank format.
            
            Comment type specific parsers (even partial) should be included
            
            comment types:
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
                ann_key='_'.join((element.tag,element.attrib['type'].replace(' ','')))
                for text_element in element.getiterator('text'):
                    if text_element.text:
                        append_to_annotations(ann_key,text_element.text)
            elif element.attrib['type']=='subcellular location':
                for subloc_element in element.getiterator('subcellularLocation'):
                    for el in subloc_element.getchildren():
                        if el.text:
                            ann_key='_'.join((element.tag,element.attrib['type'].replace(' ',''),el.tag))
                            append_to_annotations(ann_key,el.text)
            elif element.attrib['type']=='interaction':
                for interact_element in element.getiterator('interactant'):
                    ann_key='_'.join((element.tag,element.attrib['type'],'intactId'))
                    append_to_annotations(ann_key,interact_element.attrib['intactId'])
            elif element.attrib['type']=='alternative products':
                for alt_element in element.getiterator('isoform'):
                    ann_key='_'.join((element.tag,element.attrib['type'].replace(' ',''),'isoform'))
                    for id_element in alt_element.getiterator('id'):
                        append_to_annotations(ann_key,id_element.text)
            elif element.attrib['type']=='mass spectrometry':
                ann_key='_'.join((element.tag,element.attrib['type'].replace(' ','')))
                for loc_element in element.getiterator('location'):
                    pos_els=loc_element.getiterator('position')
                    if pos_els:
                        start=int(pos_els[0].attrib['position'])
                        end=start
                    else:
                        try:
                            start=int(loc_element.getiterator('begin')[0].attrib['position'])
                            end=int(loc_element.getiterator('end')[0].attrib['position'])
                        except KeyError:#undefined positions
                            start=end=0    
                mass=element.attrib['mass']
                method=element.attrib['mass']
                if start==end==0:  
                    append_to_annotations(ann_key,'undefined:%s|%s'%(mass,method))
                else:
                    append_to_annotations(ann_key,'%s..%s:%s|%s'%(start,end,mass,method))
            elif element.attrib['type']=='sequence caution':
                pass#not parsed: few information, complex structure
            elif element.attrib['type']=='online information':
                for link_element in element.getiterator('link'):
                    ann_key='_'.join((element.tag,element.attrib['type'].replace(' ','')))
                    for id_element in link_element.getiterator('link'):
                        append_to_annotations(ann_key,'%s@%s'%(element.attrib['name'],link_element.attrib['uri']))            
            
            '''return raw XML comments if needed '''
            if self.return_raw_comments:
                ann_key='_'.join((element.tag,element.attrib['type'].replace(' ',''),'xml'))
                append_to_annotations(ann_key,ElementTree.tostring(element))
                
        
        def _parse_dbReference(element):
            self.ParsedSeqRecord.dbxrefs.append(element.attrib['type']+':'+element.attrib['id'])
            for ref_element in  element.getchildren():  
                if ref_element.tag=='property':
                    pass# this data cannot be fitted in a seqrecord object with a simple list. however at least ensembl and EMBL parsing can be improved to add entries in dbxrefs
        
        def _parse_reference(element):
            reference=SeqFeature.Reference()
            authors=[]
            scopes=[]
            tissues=[]
            for ref_element in element.getchildren():
                if ref_element.tag=='citation':
                    pub_type=ref_element.attrib['type']
                    if ref_element.attrib.has_key('name'):
                        journal_name=ref_element.attrib['name']
                    else:
                        journal_name=''
                    if ref_element.attrib.has_key('date'):
                        pub_date=ref_element.attrib['date']
                    else:
                        pub_date=''
                    for cit_element in ref_element.getchildren():
                        if cit_element.tag=='title':
                            reference.title=cit_element.text
                        elif cit_element.tag=='authorList':
                            for person_element in cit_element.getchildren():
                                authors.append(person_element.attrib['name'])
                        elif cit_element.tag=='dbReference':
                            self.ParsedSeqRecord.dbxrefs.append(cit_element.attrib['type']+':'+cit_element.attrib['id'])
                            if cit_element.attrib['type']=='PubMed':
                                reference.pubmed_id=cit_element.attrib['id']
                            if ref_element.attrib['type']=='MEDLINE':
                                reference.medline_id=cit_element.attrib['id']
                elif ref_element.tag=='scope':
                    scopes.append(ref_element.text)
                elif ref_element.tag=='tissue':
                    tissues.append(ref_element.text)
            if scopes:
                scopes_str='Scopes: '+', '.join(scopes)
            else:
                scopes_str=''
            if tissues:
                tissues_str='Tissues: '+', '.join(tissues)
            else:
                tissues_str=''
            
            reference.location = [] #locations cannot be parsed since they are actually written in free text inside scopes so all the references are put in the annotation.
            reference.authors = ', '.join(authors) 
            reference.journal = journal_name 
            reference.comment = ' | '.join((pub_type,pub_date,scopes_str,tissues_str))
            append_to_annotations('references',reference)
            
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
                if feature_element.tag=='location':
                    try:
                        position_elements=feature_element.findall('position')
                        if position_elements:
                            position=int(position_elements[0].attrib['position'])
                            feature.location=SeqFeature.FeatureLocation(position,position)
                        else:
                            start_positions_elements=feature_element.findall('begin')
                            end_positions_elements=feature_element.findall('end')
                            if start_positions_elements and end_positions_elements:
                                start_position=int(start_positions_elements[0].attrib['position'])
                                end_position=int(end_positions_elements[0].attrib['position'])
                                feature.location=SeqFeature.FeatureLocation(start_position,end_position)
                    except:
                        pass#skip if parsing error    
                else:
                    try:
                        feature.qualifiers[feature_element.tag]=feature_element.text
                    except:
                        pass#skip unparsable tag
            self.ParsedSeqRecord.features.append(feature)
            
        def _parse_proteinExistence(element):
            append_to_annotations(element.tag,element.attrib['type'])   
            
        def _parse_evidence(element):
            for k, v in  element.attrib.items():
                ann_key='_'.join((element.tag,k))
                append_to_annotations(ann_key,v)   
        
        def  _parse_sequence(element):
            for k,v in element.attrib.items():
                self.ParsedSeqRecord.annotations[element.tag+'_'+k]=v
            seq=''.join((element.text.split()))
            self.ParsedSeqRecord.seq=Seq.Seq(seq,self.alphabet)
            
        #============================================#
        '''Entry attribs parsing '''
        if self.entry.attrib.has_key('dataset'):
            self.dbname=self.entry.attrib['dataset']
        else:
            self.dbname='UnknownDataset'#this should not happen!
        '''add attribs to annotations '''
        for k,v in self.entry.attrib.items():
            ann_key='_'.join((self.entry.tag,k))
            self.ParsedSeqRecord.annotations[ann_key]=v

        '''Top-to-bottom entry children parsing '''
        for element in self.entry.getchildren():
            if element.tag=='name':
                _parse_name(element)
            elif element.tag=='accession':
                _parse_accession(element)
            elif element.tag=='protein':
                _parse_protein(element)  
            elif element.tag=='gene':
                _parse_gene(element)
            elif element.tag=='geneLocation':
                _parse_geneLocation(element)
            elif element.tag=='organism':
                _parse_organism(element)          
            elif element.tag=='organismHost':    
                _parse_organismHost(element)
            elif element.tag=='keyword':
                _parse_keyword(element)
            elif element.tag=='comment':
                _parse_comment(element)
            elif element.tag=='dbReference':
                _parse_dbReference(element)
            elif element.tag=='reference':
                _parse_reference(element)
            elif element.tag=='feature':
                _parse_feature(element)
            elif element.tag=='proteinExistence':
                _parse_proteinExistence(element)
            elif element.tag=='evidence':
                _parse_evidence(element)
            elif element.tag=='sequence':
                _parse_sequence(element)
            else:
                pass   
            
        self.ParsedSeqRecord.dbxref=list(set(self.ParsedSeqRecord.dbxrefs))#remove duplicate dbxref
        self.ParsedSeqRecord.dbxref.sort()
        
        return self.ParsedSeqRecord
        