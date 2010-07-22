'''
Created on 23/apr/2010

@author: Andrea Pierleoni

'''



'''
TEST CODE
=========


::

    >>> from Bio.DAS import DASregistry 
    >>> das = DASregistry() # fetch alla available servers form DASregistry
    >>> server = das.das_servers['transmem_pred'] # choose a server
    >>> params = dict(segment = 'P50225') # set request parameters in a dictionary
    >>> results = server.fetch('das1:features', params) # fetch results
    >>> print results[0].id
    P50225
    >>> print results[0].features[0]
    type: extramembrane_region
    location: [1:295]
    id: P50225_CHAIN_1_295_CINTHIA
    strand: None
    qualifiers: 
        Key: category, Value: inferred from electronic annotation(ECO:00000067)
        Key: id, Value: SO:0001072
        Key: label, Value: P50225_CHAIN_1_295_CINTHIA
        Key: link, Value: http://pongo.biocomp.unibo.it:80/repcrc?uniprotcode=P50225
        Key: method, Value: CINTHIA
        Key: method_id, Value: CINTHIA
        
        
    #TEST a service
    >>> server._test_capab(server.capabs.keys()[0]) # uses default id to test a service
    requested  url:
    http://pongo.biocomp.unibo.it/das/dasdb/features
    
    passed params:
    segment=P50225
    
    http return code:
    200
    
    das return code:
    200
    
    http headers:
    Date: Tue, 20 Jul 2010 12:10:34 GMT
    Server: Apache/2.2.3 (Debian)
    %i 2010-07-20T08: 10:35 /das/dasdb/features 200 200
    X-DAS-Capabilities: types/2.0; features/1.0
    X-DAS-Server: ProServer/473.
    X-DAS-Status: 200
    X-DAS-Version: DAS/1.53E
    Content-Length: 665
    Last-Modified: Thu, 01 Jan 1970 00:00:00 GMT
    Connection: close
    Content-Type: text/xml
    
    
    Response: 
    
    <?xml version="1.0" standalone="yes"?>
    <!DOCTYPE DASGFF SYSTEM "http://www.biodas.org/dtd/dasgff.dtd">
    <DASGFF>
      <GFF version="1.01" href="http://localhost:9000/das/dasdb/features">
    <SEGMENT id="P50225" version="1.0" start="1" stop=""><FEATURE id="P50225_CHAIN_1_295_CINTHIA" label="P50225_CHAIN_1_295_CINTHIA"><TYPE id="SO:0001072" category="inferred from electronic annotation(ECO:00000067)">extramembrane_region</TYPE><START>1</START><END>295</END><METHOD id="CINTHIA">CINTHIA</METHOD><LINK href="http://pongo.biocomp.unibo.it:80/repcrc?uniprotcode=P50225">http://pongo.biocomp.unibo.it:80/repcrc?uniprotcode=P50225</LINK></FEATURE></SEGMENT>  </GFF>
    </DASGFF>

    
    
    
    #retrieve all available features given a uniprot ID and build a seqrecord
    
    >>> from Bio.DAS import DASregistry
    >>> das = DASregistry()
    >>> seqrec = das.fetch_to_seqrec(id = 'P00280',  
    ...                             coord_server = das.das_servers['uniprot'], 
    ...                             feat_servers = das.coords['UniProt,Protein Sequence'][:20]) # limiting to first 20 DAS servers

    >>> print seqrec
    ID: P00280
    Name: <unknown name>
    Description: <unknown description>
    Number of features: 17
    /family_annotation=['6683.1.1.1.1.1.2.1.1.1.1|Gene3D Protein Family Code: To view homologues please visit Gene3D']
    /protein secretion=['SECRETION:0:0']
    /moltype=Protein
    /version=8f29b1d09f48ec8e9909593052dcfaec
    Seq('MLAKATLAIVLSAASLPVLAAQCEATIESNDAMQYNLKEMVVDKSCKQFTVHLK...LSN', ProteinAlphabet())



    #fetch the sequence form a DAS server and return a Seq object
    >>> from Bio.DAS import DASregistry
    >>> das = DASregistry()
    >>> server = das.das_servers['uniprot']
    >>> print server.capabs.keys()
    ['das1:types', 'das1:entry_points', 'das1:features', 'das1:stylesheet', 'das1:sequence']
    >>> results = server.fetch('das1:sequence', dict(segment = 'P00280'))
    >>> dasseq = results[0]
    >>> dasseq.seq
    Seq('MLAKATLAIVLSAASLPVLAAQCEATIESNDAMQYNLKEMVVDKSCKQFTVHLK...LSN', ProteinAlphabet())


    #fetch the features form a DAS server
from Bio.DAS import DASregistry    
das = DASregistry()
server = das.das_servers['uniprot']
results = server.fetch('das1:features', dict(segment = 'P00280'))
featobj = results[0]
for feature in featobj.features[:10]:
    print feature
    
for feature in featobj.non_positional[:10]:
    print feature
    
for key,value in featobj.annotations.items():
    print key,'\t',value

'''

import urllib, urllib2
from Bio import Seq, SeqFeature, Alphabet, SeqRecord
import warnings
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
try:
    from xml.etree import cElementTree as ElementTree
except ImportError:
    try:
        from xml.etree import ElementTree as ElementTree
    except ImportError:
        # Python 2.4 -- check for 3rd-party implementations
        try:
            from lxml.etree import ElementTree
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


class DASGFFSegment(object):
    '''
    handles segments objects returned as DASGFF by DAS servers 
    '''
    def __init__(self,
                 id,
                 start,
                 stop,
                 version,
                 **kwargs):
        self.id = id
        self.start = start
        self.stop = stop
        self.version = version
        self.label = kwargs.get('label',None)
        self.type = kwargs.get('type',None)
        self.features = []
        self.non_positional = []
        self.annotations = {}
        
class DASSeq(object):
    '''
    handles sequences objects returned by DAS servers 
    '''
    def __init__(self,
                 id,
                 start,
                 stop,
                 version,
                 moltype,
                 seq):
        
        
        self.id = id
        self.start = start
        self.stop = stop
        self.version = version
        self.moltype = moltype
        seq = seq.strip()
        if self.moltype == 'DNA':
            self.seq = Seq.Seq(seq, alphabet = Alphabet.DNAAlphabet()) 
        elif  self.moltype == 'ssRNA':
            self.seq = Seq.Seq(seq, alphabet = Alphabet.RNAAlphabet()) 
        elif  self.moltype == 'dsRNA':
            self.seq = Seq.Seq(seq, alphabet = Alphabet.RNAAlphabet()) 
        elif  self.moltype == 'Protein':
            self.seq = Seq.Seq(seq, alphabet = Alphabet.ProteinAlphabet()) 
        else:
            warnings.warn('Unknown moltype: %s'%self.moltype)
            self.seq = Seq.Seq(seq,) 
            



     

class DASrequest(object):
    '''Base object to handle DAS requests'''
    
    def __init__(self, url, params = dict(), timeout = 30):
        
        self.params = urllib.urlencode(params)
        self.url = url 
        self.timeout = timeout
        self._make_call()
        
    
    def _make_call(self):
        
        self.das_response = 'unknown'
        try:
            data = urllib2.urlopen(url = self.url, 
                                   data = self.params,
                                   timeout = 30)
            if 'x-das-status' in data.headers:
                das_response = data.headers['x-das-status']
                if das_response == '402':
                    raise ValueError('server did not recognize the request')
        except ValueError:
            data = urllib2.urlopen(self.url+'?'+ self.params,
                                   timeout = 30)
            if 'x-das-status' in data.headers:
                das_response = data.headers['x-das-status']
                


        self.http_return_code = data.code
        self.das_response = das_response
        self.http_headers = data.headers
        self.response = data.read()
        
        if das_response != '200':
            if das_response == '400':
                raise Exception('DASError - Bad command (command not recognized): %s'%self.response)
            elif das_response == '401':
                raise Exception('DASError - Bad data source (data source unknown): %s'%self.response)
            elif das_response == '402':
                raise Exception('DASError - Bad command arguments (arguments invalid): %s'%self.response)
            elif das_response == '403':
                raise Exception('DASError - Bad reference object (reference sequence unknown): %s'%self.response)
            elif das_response == '404':
                raise Exception('DASError - Bad stylesheet (requested stylesheet unknown): %s'%self.response)
            elif das_response == '405':
                raise Exception('DASError - Coordinate error (sequence coordinate is out of bounds/invalid): %s'%self.response)
            elif das_response == '500':
                raise Exception('DASError - Server error, not otherwise specified: %s'%self.response)
            elif das_response == '501':
                raise Exception('DASError - Unimplemented feature: %s'%self.response)
            raise Exception('DASError - %s : %s'%(das_response, self.response))
        

class DASserver(object):
    '''
    Baseclass to represent DAS servers in the DASregistry
    docs for das at http://www.biodas.org/documents/spec.html
    '''

    def __init__(self, 
                 title,
                 description ,
                 uri ,
                 mantainer,
                 version,
                 coords,
                 capabs,
                 props,
                 elem):
        
        self.title = title
        self.description = description
        self.uri = uri
        self.mantainer = mantainer
        self.version = version
        self.coords = coords
        self.capabs = capabs
        self.props = props
        self.elem = elem


    def _make_call(self, url, params):
        return DASrequest(url, params)

    
    def _test_capab(self, capab, test=  None):
        '''test service'''
        if not test:
            try:
                coord = self.coords.keys()[0]
                test = self.coords[coord]['test_range']
            except IndexError:
                print 'no default test available'
        
        params = dict(segment = test)
        url = self.capabs[capab]['query_uri']
        response = self._make_call (url, params)

        print 'requested  url:\n', response.url
        print '\npassed params:\n', response.params
        print '\nhttp return code:\n', response.http_return_code
        print '\ndas return code:\n', response.das_response
        print '\nhttp headers:\n', response.http_headers
        print '\nResponse: \n\n', response.response
        return 
        
        
    def _parse_dasgff(self, data):
        '''parse dasgff format in biopython SeqFeatures
        '''
        
        def build_annotations(segment):
            '''try to build annotations given a list of 
            non positional features
            
            assumptions:
            annotation type == feature type
            annotation value == feature qualifier "label"  and "note" '''
            
            for feature in segment.non_positional:
                anns = []
                if 'label' in feature.qualifiers:
                    anns.append(feature.qualifiers['label'])
                if 'note' in feature.qualifiers:
                    anns.append(feature.qualifiers['note'])
                if anns:
                    if feature.type not in segment.annotations:
                        segment.annotations[feature.type] = []
                    segment.annotations[feature.type].append('|'.join(anns))
                

        element = ElementTree.parse(StringIO(data)) 
        root = element.getroot()
        
        segments = []
        if root.tag =='DASGFF':
            for seg_elem in root.getiterator('SEGMENT'):
                segment = DASGFFSegment(**seg_elem.attrib)
                for feat_elem in seg_elem.getiterator('FEATURE'):
                    feature = SeqFeature.SeqFeature(id = feat_elem.attrib['id'],
                                                    qualifiers = feat_elem.attrib)
                    for type_elem in feat_elem.getiterator('TYPE'):
                        feature.type = type_elem.text.strip()
                        feature.qualifiers.update(type_elem.attrib)
                    for start_elem in feat_elem.getiterator('START'):
                        start = int(start_elem.text)
                    for end_elem in feat_elem.getiterator('END'):
                        end = int(end_elem.text)
                    feature.location=SeqFeature.FeatureLocation(start,end)
                    for met_elem in feat_elem.getiterator('METHOD'):
                        method = met_elem.text.strip()
                        if method:
                            feature.qualifiers['method'] = method
                        if 'id' in met_elem.attrib:
                            feature.qualifiers['method_id'] = met_elem.attrib['id']
                    for link_elem in feat_elem.getiterator('LINK'):
                        feature.qualifiers['link'] = link_elem.attrib['href']
                    for score_elem in feat_elem.getiterator('SCORE'):
                        try:
                            score = float(score_elem.text.strip())
                            feature.qualifiers['score'] = score
                        except:#if '-' in score ore not valid score is provided
                            pass
                    for ori_elem in feat_elem.getiterator('ORIENTATION'):
                        ori = ori_elem.text.strip()
                        if ori == '+':
                            feature.strand = 1
                        elif ori == '-':
                            feature.strand = -1
                    for phase_elem in feat_elem.getiterator('PHASE'):
                        feature.qualifiers['phase'] = phase_elem.text.strip()
                    for note_elem in feat_elem.getiterator('NOTE'):
                        feature.qualifiers['note'] = note_elem.text.strip()  
                    for target_elem in feat_elem.getiterator('TARGET'):
                        feature.qualifiers['target_id'] = target_elem.attrib['id']
                        feature.qualifiers['target_start'] = target_elem.attrib['start']
                        feature.qualifiers['target_stop'] = target_elem.attrib['stop']
                    for group_elem in feat_elem.getiterator('GROUP'):
                        feature.qualifiers['group_id'] = group_elem.attrib['id']
                        if 'label' in group_elem.attrib:
                            feature.qualifiers['group_label'] = group_elem.attrib['label']
                        if 'type' in group_elem.attrib:
                            feature.qualifiers['group_type'] = group_elem.attrib['type']
                    
                    if start == end == 0: #non positional feature
                        segment.non_positional.append(feature)
                    else:    
                        segment.features.append(feature)
                
                build_annotations(segment)
                segments.append(segment)
        else:
            raise ValueError('data not in valid DASGFF format')

        return segments
    
    def _parse_dassequence(self, data):
        '''parse DASSEQUENCE responses in biopython Seq objects
        <?xml version="1.0" standalone='no'?>
        <!DOCTYPE DASSEQUENCE SYSTEM "http://www.biodas.org/dtd/dassequence.dtd">
        <DASSEQUENCE>
          <SEQUENCE id="P00280" start="1" stop="149" moltype="Protein" version="8f29b1d09f48ec8e9909593052dcfaec">MLAKATLAIVLSAASLPVLAAQCEATIESNDAMQYNLKEMVVDKSCKQFTVHLKHVGKMAKVAMGHNWVLTKEADKQGVATDGMNAGLAQDYVKAGDTRVIAHTKVIGGGESDSVTFDVSKLTPGEAYAYFCSFPGHWAMMKGTLKLSN</SEQUENCE>
        </DASSEQUENCE>
        '''
        
        element = ElementTree.parse(StringIO(data)) 
        root = element.getroot()
        
        sequences = []
        if root.tag =='DASSEQUENCE':
            for seq_elem in root.getiterator('SEQUENCE'):
                sequences.append(DASSeq(id = seq_elem.attrib['id'],
                                         start= seq_elem.attrib['start'],
                                         stop= seq_elem.attrib['stop'],
                                         version = seq_elem.attrib['version'],
                                         moltype = seq_elem.attrib['moltype'],
                                         seq = seq_elem.text))
                
        return sequences

    def _parse_dasstyle(self, data):
        '''<DASSTYLE>  <STYLESHEET>    <CATEGORY id="default">      <TYPE id="default">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#666666</FGCOLOR>                <BGCOLOR>#CCCCCC</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="ACT_SITE">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#000000</FGCOLOR>                <BGCOLOR>#111111</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="BINDING">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#584858</FGCOLOR>                <BGCOLOR>#8B7B8B</BGCOLOR>        </BOX>        </GLYPH>     </TYPE>     <TYPE id="CA_BIND">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CCCC00</FGCOLOR>                <BGCOLOR>#ffff00</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="CARBOHYD">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#581456</FGCOLOR>                <BGCOLOR>#8B4789</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="CHAIN">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC7200</FGCOLOR>                <BGCOLOR>#FFA500</BGCOLOR>            </BOX>    </GLYPH>     </TYPE>     <TYPE id="component">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#0000CC</FGCOLOR>                <BGCOLOR>#3333ff</BGCOLOR>            </BOX>    </GLYPH>     </TYPE>     <TYPE id="CONFLICT">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>    <FGCOLOR>#1A1A1A</FGCOLOR>                <BGCOLOR>#4D4D4D</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="CROSSLINK">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#9A520C</FGCOLOR>                <BGCOLOR>#CD853F</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="DISULFID">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#00CC00</FGCOLOR>                <BGCOLOR>#00FF00</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="DNA_BIND">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#58002F</FGCOLOR>                <BGCOLOR>#8B1C62</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="DOMAIN">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#9A3300</FGCOLOR>                <BGCOLOR>#CD6600</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="HELIX">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#005800</FGCOLOR>                <BGCOLOR>#008B00</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="INIT_MET">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#0000CC</FGCOLOR>                <BGCOLOR>#0000FF</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="LIPID">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BBBB00</FGCOLOR>                <BGCOLOR>#EEEE00</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="METAL">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#033158</FGCOLOR>                <BGCOLOR>#36648B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="MOD_RES">        <GLYPH        <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#37279A</FGCOLOR>                <BGCOLOR>#6A5ACD</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="MUTAGEN">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC8D98</FGCOLOR>                <BGCOLOR>#FFC0CB</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="NON_CONS">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#003100</FGCOLOR>                <BGCOLOR>#006400</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="NON_TER">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#000058</FGCOLOR>                <BGCOLOR>#00008B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="NP_BIND">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#580300</FGCOLOR>                <BGCOLOR>#8B3626</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="PEPTIDE">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#475858</FGCOLOR>                <BGCOLOR>#7A8B8B</BGCOLOR>        </BOX>        </GLYPH>     </TYPE>     <TYPE id="Pfam">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#3300CC</FGCOLOR>                <BGCOLOR>#6633ff</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="PIR">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC0000</FGCOLOR>                <BGCOLOR>#ff3300</BGCOLOR>            </BOX>    </GLYPH>     </TYPE>     <TYPE id="PRINTS">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>    <FGCOLOR>#009900</FGCOLOR>                <BGCOLOR>#33cc33</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="ProDom">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#3366CC</FGCOLOR>                <BGCOLOR>#6699ff</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="PROPEP">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB3774</FGCOLOR>                <BGCOLOR>#EE6AA7</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="PROSITE">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC6600</FGCOLOR>                <BGCOLOR>#ff9933</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="REPEAT">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#580000</FGCOLOR>                <BGCOLOR>#8B1A1A</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="Scop">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#666600</FGCOLOR>                <BGCOLOR>#999900</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SE_CYS">        <GLYPH>        <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC00CC</FGCOLOR>                <BGCOLOR>#FF00FF</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SIGNAL">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB7BBB</FGCOLOR>                <BGCOLOR>#EEAEEE</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SIMILAR">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#36269A</FGCOLOR>                <BGCOLOR>#6959CD</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SITE">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#005858</FGCOLOR>                <BGCOLOR>#008B8B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="Smart">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#9B0000</FGCOLOR>                <BGCOLOR>#ce0031</BGCOLOR>        </BOX>        </GLYPH>     </TYPE>     <TYPE id="STRAND">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB8181</FGCOLOR>                <BGCOLOR>#EEB4B4</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="THIOETH">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#000066</FGCOLOR>                <BGCOLOR>#000099</BGCOLOR>            </BOX>    </GLYPH>     </TYPE>     <TYPE id="THIOLEST">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>    <FGCOLOR>#000066</FGCOLOR>                <BGCOLOR>#333399</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="TIGRFAMs">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#330000</FGCOLOR>                <BGCOLOR>#660033</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="TRANSIT">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#9A520C</FGCOLOR>                <BGCOLOR>#CD853F</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="TRANSMEM">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#365836</FGCOLOR>                <BGCOLOR>#698B69</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="TURN">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#5E00BB</FGCOLOR>                <BGCOLOR>#912CEE</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="UNSURE">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#333333</FGCOLOR>                <BGCOLOR>#000000</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="VARIANT">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#330000</FGCOLOR>                <BGCOLOR>#660033</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="VARSPLIC">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB4FBB</FGCOLOR>                <BGCOLOR>#EE82EE</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="ZN_FING">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC00CC</FGCOLOR>                <BGCOLOR>#FF00FF</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SMART">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#9B0000</FGCOLOR>                <BGCOLOR>#ce0031</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="ISOFORM">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB4FBB</FGCOLOR>                <BGCOLOR>#EE82EE</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001104">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#000000</FGCOLOR>                <BGCOLOR>#111111</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001091">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#584858</FGCOLOR>                <BGCOLOR>#8B7B8B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0000417">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CCCC00</FGCOLOR>                <BGCOLOR>#ffff00</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="MOD:00693">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#581456</FGCOLOR>                <BGCOLOR>#8B4789</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0000419">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC7200</FGCOLOR>                <BGCOLOR>#FFA500</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001085">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#1A1A1A</FGCOLOR>                <BGCOLOR>#4D4D4D</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001087">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#9A520C</FGCOLOR>                <BGCOLOR>#CD853F</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001088">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#00CC00</FGCOLOR>                <BGCOLOR>#00FF00</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001117">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#005800</FGCOLOR>                <BGCOLOR>#008B00</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0000691">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#0000CC</FGCOLOR>                <BGCOLOR>#0000FF</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001092">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#033158</FGCOLOR>                <BGCOLOR>#36648B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001089">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#37279A</FGCOLOR>                <BGCOLOR>#6A5ACD</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001148">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC8D98</FGCOLOR>                <BGCOLOR>#FFC0CB</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001083">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#003100</FGCOLOR>                <BGCOLOR>#006400</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001084">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#000058</FGCOLOR>                <BGCOLOR>#00008B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001064">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#475858</FGCOLOR>                <BGCOLOR>#7A8B8B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001062">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB3774</FGCOLOR>                <BGCOLOR>#EE6AA7</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001068">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#580000</FGCOLOR>                <BGCOLOR>#8B1A1A</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="MOD:00031">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#CC00CC</FGCOLOR>                <BGCOLOR>#FF00FF</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0000418">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB7BBB</FGCOLOR>                <BGCOLOR>#EEAEEE</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0000839">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#005858</FGCOLOR>                <BGCOLOR>#008B8B</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001111">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#BB8181</FGCOLOR>                <BGCOLOR>#EEB4B4</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0000725">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#9A520C</FGCOLOR>                <BGCOLOR>#CD853F</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001077">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#365836</FGCOLOR>                <BGCOLOR>#698B69</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001128">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#5E00BB</FGCOLOR>                <BGCOLOR>#912CEE</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001086">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#333333</FGCOLOR>                <BGCOLOR>#000000</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>     <TYPE id="SO:0001147">        <GLYPH>            <BOX>                <LINEWIDTH>1</LINEWIDTH>                <FGCOLOR>#330000</FGCOLOR>                <BGCOLOR>#660033</BGCOLOR>            </BOX>        </GLYPH>     </TYPE>   </CATEGORY>  </STYLESHEET></DASSTYLE>'''
        return 'TO DO'
    
    def _parse_dastypes(self, data):
        '''<?xml version="1.0" standalone='no'?>
            <!DOCTYPE DASTYPES SYSTEM "http://www.biodas.org/dtd/dastypes.dtd">
            <DASTYPES>
              <GFF version="1.0" href="http://www.ebi.ac.uk/das-srv/uniprot/das/uniprot/types">
                <SEGMENT version="2.0" label="Complete datasource summary">
                  <TYPE id="lipoconjugated residue" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="sequence_caution" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="subunit_information" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="sequence_conflict" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="disease" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="compositionally_biased_region" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="annotation_warnings" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="biotechnological_use" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="non_adjacent_residue" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="L-selenocysteine" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="crosslinked residues" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="cofactor" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="metal_contact" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="non_terminal_residue" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="binding_motif" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="polypeptide_repeat" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="cleaved_initiator_methionine" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="alpha_helix" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="mutated_variant_site" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="taxonomy" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="coiled_coil" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="web_resource" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="toxic_dose" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="biological_process" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="turn" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="polypeptide_motif" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="cellular_component" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="pathway" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="alternative_products" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="EC_annotation" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="catalytic_activity" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="glycosylated residue" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="tissue_specificity" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="extramembrane" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="sequence_uncertainty" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="post_translational_modification" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="protein_name" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="allergen" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="induction" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="developmental_stage" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="post_translational_modification_information" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="biophysiochemical_properties" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="organism_species" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="enzyme_regulation" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="polypeptide_region" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="beta_strand" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="publication" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="mass_spectrometry" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="active_peptide" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="propeptide" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="natural_variant_site" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="catalytic_residue" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="polymorphism" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="polypeptide_domain" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="technical_term" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="ligand_information" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="transit_peptide" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="mature_protein_region" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="alternative_sequence_site" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="miscellaneous" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="pharmaceutical" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="functional_annotation" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="RNA_editing" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="protein_protein_interactions" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="signal_peptide" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="transmembrane" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="disulfide crosslinked residues" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="domain_information" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="similarity" method="UniProt" category="inferred by curator (ECO:0000001)" />
                  <TYPE id="molecular_function" method="UniProt" category="inferred by curator (ECO:0000001)" />
                </SEGMENT>
              </GFF>
            </DASTYPES>
            '''                     
        return 'TO DO'
    
    def _parse_dasep(self, data):
        '''
        <?xml version="1.0" standalone='no'?>
        <!DOCTYPE DASEP SYSTEM "http://www.biodas.org/dtd/dasep.dtd">
        <DASEP>
          <ENTRY_POINTS href="http://www.ebi.ac.uk/das-srv/uniprot/das/uniprot/entry_points" version="2010.05">
            <SEGMENT id="A0" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A1" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A2" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A3" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A4" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A5" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A6" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A7" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A8" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="A9" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B0" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B1" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B2" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B3" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B4" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B5" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B6" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B7" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B8" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="B9" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C0" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C1" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C2" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C3" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C4" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C5" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C6" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C7" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C8" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="C9" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="D0" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="D1" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="D2" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="D3" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="D4" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O0" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O1" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O2" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O3" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O4" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O5" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O6" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O7" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O8" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="O9" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P0" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P1" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P2" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P3" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P4" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P5" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P6" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P7" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P8" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="P9" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q0" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q1" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q2" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q3" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q4" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q5" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q6" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q7" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q8" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
            <SEGMENT id="Q9" start="1" stop="1" orientation="+" subparts="yes">Scaffold Entry Point</SEGMENT>
          </ENTRY_POINTS>
        </DASEP>

        '''
        return 'TO DO'
    
    def fetch(self, capab, params,):
        url = self.capabs[capab]['query_uri']
        response = self._make_call (url, params)
        data = response.response
        data_root = ElementTree.parse(StringIO(data)).getroot().tag
        if data_root == 'DASGFF':
            return self._parse_dasgff(data)
        elif data_root == 'DASSEQUENCE':
            return self._parse_dassequence(data)
        elif data_root == 'DASSTYLE':
            return self._parse_dasstyle(data)
        elif data_root == 'DASTYPES':
            return self._parse_dasgff(data)
        elif data_root == 'DASEP':
            return self._parse_dasep(data)
        else:
            raise ValueError('format not supported: %s'%data_root)
        
    
    
    def _read_dsn(self):
        '''if needed parse the list of data sources for this DAS server
        
        url to parse: 
        
        Retrieve the List of Data Sources
        
        Scope: Reference and annotation servers.
        
        Command: dsn
        
        Format:
        
            PREFIX/das/dsn
        
        Description: This query returns the list of data sources that are available from this server.
        Response:
        
        The response to the dsn command is the "DASDSN" XML-formatted document:
        
        Format:
        
            <?xml version="1.0" standalone="no"?>
            <!DOCTYPE DASDSN SYSTEM "http://www.biodas.org/dtd/dasdsn.dtd">
            <DASDSN>
              <DSN>
                <SOURCE id="id1" version="version">source name 1</SOURCE>
                <MAPMASTER>URL</MAPMASTER>
                <DESCRIPTION>descriptive text 1</DESCRIPTION>
              </DSN>
              <DSN>
                <SOURCE id="id2" version="version">source name 2</SOURCE>
                <MAPMASTER>URL</MAPMASTER>
                <DESCRIPTION href="url">descriptive text 2</DESCRIPTION>
              </DSN>
              ...
            </DASDSN>
        
        <!DOCTYPE> (required; one only)
            The doctype indicates which formal DTD specification to use. For the dsn query, the doctype DTD is "http://www.biodas.org/dtd/dasdsn.dtd".
        
        <DASDSN> (required; one only)
            The appropriate doctype and root tag is DASDSN.
        
        <DSN> (required; one or more)
            There are one or more <DSN> tags, one for each data source. Each <DSN> contains one <SOURCE> tag, one <MAPMASTER> tag, and optionally one <DESCRIPTION> tag.
        
        <SOURCE> (required; one per DSN tag)
            This tag indicates the symbolic name for a data source. The symbolic name to use for further requests can be found in the id (required) attribute. A source version attribute is optional, but strongly recommended. The tag body contains a human-readable label which may or may not be different from the ID.
        
        <MAPMASTER> (required; one per DSN tag)
            This tag contains the URL (site.specific.prefix/das/data_src) that is being annotated by this data source. For an annotation server, this is the reference server which is being annotated. For a reference server, this would echo its own URL.
        
        <DESCRIPTION> (optional)
            This tag contains additional descriptive information about the data source. If an href (optional) attribute is present, the attribute contains a link to further human-readable information about the data source, such as its home page. 
        '''
        return 'TO DO'
    



class DASregistry(object):
    '''
    Reads sources from DASregistry and returns links to DAS servers
    '''


    def __init__(self, ):
        '''
        fetch all the available sources
        and returns a connection object 
        
        >>> das = DASregistry()
        
        with this properties:
        
        - das.auths       ---> DAS Authorities
        - das.coords      ---> DASserver titles grouped by coordinate system
        - das.das_servers ---> dictionary containing alla available DASserver objects
        '''
        self.das_servers=dict()
        sources = urllib.urlopen('http://www.dasregistry.org/das/sources')
        for event, elem in ElementTree.iterparse(sources):
            if event=="end" and elem.tag =='SOURCE':
                try:
                    server_title = elem.attrib['title']
                    coords = {}
                    capabs = {}
                    props = []
                    for man_elem in  elem.getiterator('MAINTAINER'):
                        mantainer = man_elem.attrib['email']
                    for ver_elem in  elem.getiterator('VERSION'):
                        version = ver_elem.attrib
                    for cor_elem in  elem.getiterator('COORDINATES'):
                        cor_value  = cor_elem.text.strip()
                        coords[cor_value] = cor_elem.attrib
                    for cap_elem in  elem.getiterator('CAPABILITY'):
                        cap_value  = cap_elem.attrib['type']
                        capabs[cap_value] = cap_elem.attrib
                    for prop_elem in  elem.getiterator('PROP'):
                        props.append(prop_elem.attrib)
                    
                    
                    self.das_servers[server_title] = DASserver(title = server_title,
                                                              description = elem.attrib['description'],
                                                              uri = elem.attrib['uri'],
                                                              mantainer = mantainer,
                                                              version = version,
                                                              coords = coords,
                                                              capabs = capabs,
                                                              props = props,
                                                              elem = elem)
                except Exception, error:
                    print error
        
        #populate authorities
        self.coords = dict()#coordinates 2 server dict
        self.auths = dict()#authorities tree
        for server_title, server in self.das_servers.items():
            for coord, auth in server.coords.items():
                if coord not in self.coords:
                    self.coords[coord] =[]
                self.coords[coord].append(server_title)
                if auth['authority'] not in self.auths:
                    self.auths[auth['authority']] = []
                if auth['source'] not in self.auths[auth['authority']]:
                    self.auths[auth['authority']].append(auth['source'])
                


    def fetch_to_seqrec(self,
                        id,  # id to be searched in the coordinate server
                        coord_server, # coordinate server used to fetch sequence
                        feat_servers, #servers to fetch features from
                        start = 'not set',
                        stop = 'not set',): 
        ''' Helper function - Fetchs all the available features given coordinates
            and return a SeqRecord object'''

        def merge_annotations(ann1, ann2):
            for k in ann2:
                if k in ann1:
                    ann1[k].extend(ann2[k])
                    ann1[k]=list(set(ann1[k]))
                else:
                    ann1[k] = ann2[k]
            return ann1

        if 'das1:sequence' in coord_server.capabs:
            sequence = coord_server.fetch('das1:sequence', dict(segment = id))
            if sequence:
                sequence = sequence[0]
                seqrec = SeqRecord.SeqRecord(sequence.seq, id = id)
                seqrec.annotations['moltype'] = sequence.moltype
                seqrec.annotations['version'] = sequence.version
            else:
                raise IndexError('sequence not found in the specified DAS server')
        else:
            raise IndexError('specified DAS server do not accept das1:sequence method')
        
        params = dict(segment = id)
        for das_server_title in feat_servers:
            das_server = self.das_servers[das_server_title]
            if 'das1:features' in das_server.capabs:
                try:
                    feature_result = das_server.fetch('das1:features', params)
                except Exception,error:
                    warnings.warn('Failed searching features for "%s" DAS server. Error: %s'%(das_server_title,repr(error)))
                        
                if feature_result:
                    if isinstance(feature_result[0],DASGFFSegment):
                        featobj = feature_result[0]
                        if featobj.features:
                            seqrec.features.extend(featobj.features)
                        if featobj.annotations:
                            seqrec.annotations = merge_annotations(seqrec.annotations, featobj.annotations)
                    
                    else:
                        warnings.warn('"%s" DAS server returned a "%s" object and not a DASGFFSegment'%(das_server_title,type(feature_result[0])))
        
        return seqrec
        
