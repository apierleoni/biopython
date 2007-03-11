# Copyright 2006, 2007 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Alphabet import generic_alphabet

class SequenceIterator :
    """Base class for building Sequence iterators.

    You should write a next() method to return SeqRecord
    objects.  You may wish to redefine the __init__
    method as well.
    """
    def __init__(self, handle, alphabet=generic_alphabet) :
        """Create a SequenceIterator object.

        handle - input file
        alphabet - optional, e.g. Bio.Alphabet.generic_protein

        Note when subclassing:
        - there should be a single non-optional argument,
          the handle.
        - you do not have to require an alphabet.
        - you can add additional optional arguments."""
        self.handle = handle
        self.alphabet = alphabet
        #####################################################
        # You may want to subclass this, for example        #
        # to read through the file to find the first record,#
        # or if additional arguments are required.          #
        #####################################################

    def next(self) :
        """Return the next record in the file
        This method should be replaced by any derived class to do something useful."""
        raise NotImplementedError, "This object should be subclassed"
        #####################################################
        # You SHOULD subclass this, to split the file up    #
        # into your individual records, and convert these   #
        # into useful objects, e.g. return SeqRecord object #
        #####################################################

    def __iter__(self):
        """Iterate over the entries as a SeqRecord objects

        Example usage for Fasta files:

        myFile = open("example.fasta","r")
        myFastaReader = FastaIterator(myFile)
        for record in FastaIterator :
            print record.id
            print record.seq
        myFile.close()"""
        return iter(self.next, None)

    def close(self):
        """Close the input file handle"""
        return self.handle.close()

class InterlacedSequenceIterator(SequenceIterator) :
    """This class should be subclassed by any iterator for a non-sequential file type.

    This object is not intended for use directly.
    
    When writing a parser for any interlaced sequence file where the whole
    file must be read in order to extract any single record, then you should
    subclass this object.

    All you need to do is to define your own:
    (1) __init__ method to parse the file and call self.move_start()
    (2) __len__ method to return the number of records
    (3) __getitem__ to return any requested record.

    This class will then provide the iterator methods including next(), but relies
    on knowing the total number of records and tracking the pending record index in
    as self._n

    It is up to the subclassed object to decide if it wants to generate a cache of
    SeqRecords when initialised, or simply use its own lists and dicts and create
    SeqRecords on request.
    """

    def __init__(self) :
        """Create the object

        This method should be replaced by any derived class to do something useful."""
        #We assume that your implementation of __init__ will ensure self._n=0
        self.move_start()
        raise NotImplementedError, "This object method should be subclassed"
        #####################################################
        # You SHOULD subclass this                          #
        #####################################################

    def __len__(self) :
        """Return the number of record

        This method should be replaced by any derived class to do something useful."""
        raise NotImplementedError, "This object method should be subclassed"
        #####################################################
        # You SHOULD subclass this                          #
        #####################################################

    def __getitem__(self, i) :
        """Return the requested record

        This method should be replaced by any derived class to do something
        useful.

        It should NOT touch the value of self._n"""
        raise NotImplementedError, "This object method should be subclassed"
        #####################################################
        # You SHOULD subclass this                          #
        #####################################################

    def move_start(self) :
        self._n = 0

    def next(self) :
        next_record = self._n
        if next_record < len(self) :
            self._n = next_record+1
            return self[next_record]
        else :
            #StopIteration
            return None
    
    def __iter__(self):
        return iter(self.next, None)

class SequenceWriter:
    """This class should be subclassed.

    Interlaced file formats (e.g. Clustal) should subclass directly.

    Sequential file formats (e.g. Fasta, GenBank) should subclass
    the SequentialSequenceWriter class instead.
    """
    def __init__(self, handle):
        """Creates the writer object

        Use the method write_file() to actually record your sequence records."""
        self.handle = handle
    
    def write_file(self, records) :
        """Use this to write an entire file containing the given records.

        records - A list or iterator returning SeqRecord objects

        This method can only be called once."""
        #Note when implementing this, you should close the file at the end.
        raise NotImplementedError, "This object should be subclassed"
        #####################################################
        # You SHOULD subclass this                          #
        #####################################################

class SequentialSequenceWriter(SequenceWriter):
    """This class should be subclassed.

    It is intended for sequential file formats with an (optional)
    header, repeated records, and an (optional) footer.

    In this case (as with interlaced file formats), the user may
    simply call the write_file() method and be done.

    However, they may also call the write_header(), followed
    by mulitple calls to write_record() and/or write_records()
    followed by write_footer() and close()

    Users must call write_header() and write_footer() even when
    the file format concerned doesn't have a header or footer.
    This is to try and make life as easy as possible when
    switching the output format.
    
    Note that write_header() cannot require any assumptions about
    the number of records - which would be required for phylip files
    for example.
    """
    def __init__(self, handle):
        self.handle = handle
        self._header_written = False
        self._record_written = False
        self._footer_written = False

    def write_header(self) :
        assert not self._header_written, "You have aleady called write_header()"
        assert not self._record_written, "You have aleady called write_record() or write_records()"
        assert not self._footer_written, "You have aleady called write_footer()"
        self._header_written = True
        
    def write_footer(self) :
        assert self._header_written, "You must call write_header() first"
        assert self._record_written, "You have not called write_record() or write_records() yet"
        assert not self._footer_written, "You have aleady called write_footer()"
        self._footer_written = True
        
    def write_record(self, record):
        """Write a single record to the output file.

        record - a SeqRecord object

        Once you have called write_header() you can call write_record()
        and/or write_records() as many times as needed.  Then call
        write_footer() and close()."""
        assert self._header_written, "You must call write_header() first"
        assert not self._footer_written, "You have already called write_footer()"
        self._record_written = True
        raise NotImplementedError, "This object should be subclassed"
        #####################################################
        # You SHOULD subclass this                          #
        #####################################################
        
    def write_records(self, records):
        """Write multiple record to the output file.

        records - A list or iterator returning SeqRecord objects

        Once you have called write_header() you can call write_record()
        and/or write_records() as many times as needed.  Then call
        write_footer() and close()."""
        #Default implementation:
        assert self._header_written, "You must call write_header() first"
        assert not self._footer_written, "You have already called write_footer()"
        for record in records :
            self.write_record(record)
        #Mark as true, even if there where no records
        self._record_written = True

    def write_file(self, records) :
        """Use this to write an entire file containing the given records.

        records - A list or iterator returning SeqRecord objects

        This method can only be called once."""
        self.write_header()
        self.write_records(records)
        self.write_footer()
        #Don't automatically close the handle:
        #self.handle.close()

    def flush(self):
        """Flush any output pending on the output file handle"""
        return self.handle.flush()

    def close(self):
        """Close the output file handle"""
        assert self._header_written, "You must call write_header() first"
        assert self._record_written, "You must call write_record() or write_records() first"
        assert self._footer_written, "You must call write_footer() first"
        return self.handle.close()
