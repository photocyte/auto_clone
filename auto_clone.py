# -*- coding: UTF-8 -*-
import pandas
import Bio
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.Alphabet
import time
import os

class IDT_plate_object():
    def __init__(self,PrimerNamePrefix,PrimerIndex,StartWellPosition,direction,ignoreWellOverlap=True):
        self.PrimerNamePrefix = PrimerNamePrefix
        self.PrimerIndex = PrimerIndex
        self.WellPosition = StartWellPosition
        self.plate_direction = direction ##alphabetical (A1,B1,C1... A2,B2,C2...) or numerical (A1,A2,A3...B1,B2,C3...)
        self.usedWells = []
        self.df = pandas.DataFrame(columns=["WellPosition","Name","Sequence","Notes"]) ##Write the header
        self.ignoreWellOverlap = ignoreWellOverlap        

    def incrementWellPositionCol(self,increment):
        self.usedWells.append(self.WellPosition)
        letter = self.WellPosition[0]
        number = int(self.WellPosition[1:])
        number += increment
        
        ##Handle wrapping
        if number > 12: ##If the number is too large
            while number > 12:
                letter = chr(ord(letter) + 1)
                number -= 12     
        if number < 1:  ##If the number is too small
            while number < 1:
                letter = chr(ord(letter) - 1)
                number += 12
        
        
        if ord(letter) > ord("I"):
            print("plate well out of range: too high")
        if ord(letter) < ord("A"):
            print("plate well out of range: too low")
        self.WellPosition = letter+str(number)
        
    def incrementWellPositionRow(self,increment):
        self.usedWells.append(self.WellPosition)
        letter = self.WellPosition[0]
        number = int(self.WellPosition[1:])
        
        ##No wrapping handling implemented.  This is because running off the letters of the plate is probably not intended (have to implement having multiple plates)
        letter = chr(ord(letter) + increment)
        if ord(letter) > ord("I"):
            print("plate well out of range: too high")
            print("Exiting program...")
            exit()
        if ord(letter) < ord("A"):
            print("plate well out of range: too low")
            print("Exiting program...")
            exit()
        self.WellPosition = letter+str(number)
        
    def incrementWellPosition(self,increment):
        if self.plate_direction == "numerical":
            self.incrementWellPositionCol(increment)
        elif self.plate_direction == "alphabetical":
            self.incrementWellPositionRow(increment)
            print("Alphabetical incrementing not properly implemented!")
            exit()
        else:
            print("Plate direction not properly set!")
            exit()
                    
    def formatPrimerRow(self,record,plasmid,direction):
        newRow = dict()
        newRow['WellPosition'] = self.WellPosition[0]+str(self.WellPosition[1:])
        if direction == "+":
            newRow['Notes'] = "Forward "+plasmid.description+record.id+record.description
            newRow['Sequence'] = plasmid.gibson_forward_overlap+str(record.Sequence)

        if direction == "-":
            newRow['Notes'] = "Reverse "+plasmid.description+record.id+record.description
            newRow['Sequence'] = plasmid.gibson_reverse_overlap+str(record.Sequence.reverse_complement())
        newRow['Name'] = Name
        return newRow
        
        
    def newPrimer(self,record):
        if self.WellPosition in self.usedWells and self.ignoreWellOverlap == False:
            print("**Primer collision** Tried to put a primer in a well that was already marked as used.")
            print("Exiting...")
            exit()
        
        if record.description == None:
            record.description = str(len(record.seq))
        
        ##Layout is ["WellPosition","Name","Sequence","Notes"]
        newRow = {'WellPosition':self.WellPosition,"Name":str(self.PrimerNamePrefix)+format(self.PrimerIndex,'04'),"Sequence":str(record.seq),"Notes":str(record.description)}
        self.df = self.df.append(newRow,ignore_index=True)
        self.usedWells = list(self.df['WellPosition'].values)
        self.PrimerIndex += 1
    
    def newPrimerPairVertical(self,record_1,record_2):
        self.newPrimer(record_1) ##Add primer A1
        self.incrementWellPositionRow(1) ##+1 for the letter A1->B1
        self.newPrimer(record_2) ##add primer B1

        self.incrementWellPositionRow(-1) ## B1->A1
    
        if int(self.WellPosition[1:]) == 12: ##Wrapping support.  Go an extra row to not run into R primer of previous primers.  E.g, if at A12, go to B12, then wrap to C1 when incremeting
            self.incrementWellPositionRow(1)
            self.incrementWellPositionCol(1) ##A1->A2 
            
    def getCSV(self):
        ##Need name, sequence, scale, purification
        ##Scale = 25nm
        ##Purification = STD
        rows = self.df.values.tolist()
        buffer = ''
        for row in rows:
            name = row[1]
            sequence = row[2]
            scale = "25nm"
            purity = "STD"
            buffer += name+","+sequence+","+scale+","+purity+'\n'

        return buffer
        
class plasmid_object():
    def __init__(self,name):
        if name != None:
            ###Load plasmids to clone into:
            ###Expects the plasmid to have a stretch of "N"s, which is where it will stick the insert.
            plasmid_path="../../plasmids/"+name+".gb"
            if not os.path.isfile(plasmid_path):
                print("No plasmid file found:",plasmid_path)
                self.plasmid = None
                self.forward_plasmid = None
                self.reverse_plasmid = None
            elif os.path.isfile(plasmid_path):
                handle = open(plasmid_path,"rU")
                self.plasmid = list(Bio.SeqIO.parse(handle,"gb"))[0]
                handle.close()
                if len(self.plasmid.seq) == 0:
                    print("Empty plasmid.")
                    self.plasmid = None
                    self.forward_plasmid = None
                    self.reverse_plasmid = None
                elif len(self.plasmid.seq) > 0:
                    self.plasmid.id = name
                    Nstart = self.plasmid.seq.find("N")
                    Nend = self.plasmid.seq.rfind("N")+1
                    self.forward_plasmid = self.plasmid[0:Nstart]
                    self.reverse_plasmid = self.plasmid[Nend:]
                
            ##Load the plasmid overlap & comments.  Also possible to derive from plasmid file.
            text_path = "../../plasmids/"+name+".txt"
            handle = open(text_path,"rU")
            forwardString = handle.readline() ##Line 1
            reverseString = handle.readline() ##Line 2
            commentString = handle.readline() ##Line 3
            handle.close()
            self.gibson_forward_overlap = forwardString[forwardString.find(":")+1:].strip()
            self.gibson_reverse_overlap = reverseString[reverseString.find(":")+1:].strip()    ##Overlap Not including stop codon
            self.description = commentString[commentString.find(":")+1:].strip()
        else:
            print("Couldn't match plasmid file!")
            self.plasmid = None
            self.forward_plasmid = None
            self.reverse_plasmid = None
            self.gibson_forward_overlap = None 
            self.gibson_reverse_overlap = None
            self.description = "uninitialized plasmid"
    
    ##This function concatenates multiple records to form a final plasmid & writes it to disk (genbank format)
    def writePlasmidFile(self,record,plate=None):
    
        if self.plasmid is None:
            print("Can't write. Plasmid wasn't loaded.")
            return None
            
        ##Have to get the record sequence type to IUPACAmbiguousDNA
        seq_reformatted = Bio.Seq.Seq(str(record.seq), Bio.Alphabet.IUPAC.IUPACAmbiguousDNA())

        record_reformatted = Bio.SeqRecord.SeqRecord(seq_reformatted,id=record.id, name=record.name,description=record.description)
        ##Have to add a feature so it looks nice in APE & other plasmid editors.  This site https://www.biostars.org/p/57549/ explains how.
        ##http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html also useful.
        my_start_pos = Bio.SeqFeature.ExactPosition(0)
        my_end_pos = Bio.SeqFeature.ExactPosition(len(record.seq))
        my_feature_location = Bio.SeqFeature.FeatureLocation(my_start_pos,my_end_pos)
        my_feature_type = "CDS"
        
        if plate != None:
            my_label = plate.WellPosition+"_"+record.id+" CDS "+str(len(record.seq))+"bp"
        else:
            my_label = record.id+" CDS "+str(len(record.seq))+"bp"
            
        my_feature = Bio.SeqFeature.SeqFeature(my_feature_location,type=my_feature_type,strand=1,qualifiers={"label":my_label,"ApEinfo_fwdcolor":"#00ccff","ApEinfo_revcolor":"#00ccff"})
        record_reformatted.features.append(my_feature)

        complete_plasmid = self.forward_plasmid+record_reformatted+self.reverse_plasmid ##Concatenated records is as simple as addition. Thanks BioPython!
        
        if plate != None:
            file_string = plate.WellPosition+"_"+record.id+"_"+self.forward_plasmid.id+".gb"
        else:
            file_string = record.id+"_"+self.forward_plasmid.id+".gb"
            
        Bio.SeqIO.write(complete_plasmid, file_string, "gb")
        
        return file_string
        
    def formatGibsonPrimerForward(self,record,forwardPrimerLength):
        primerDescription = "Forward "+self.description+record.id+record.description
        sequence = str(self.gibson_forward_overlap)+str(record.seq[:forwardPrimerLength])
        
        seq_reformatted = Bio.Seq.Seq(str(sequence), Bio.Alphabet.IUPAC.IUPACAmbiguousDNA())
        newRecord = Bio.SeqRecord.SeqRecord(seq_reformatted,id=record.id, name=record.name,description=primerDescription)
        return newRecord

    def formatGibsonPrimerReverse(self,record,reversePrimerLength):
        primerDescription = "Reverse "+self.description+record.id+record.description
        sequence = str(self.gibson_reverse_overlap)+str(record.seq[-reversePrimerLength:].reverse_complement())
        
        seq_reformatted = Bio.Seq.Seq(str(sequence), Bio.Alphabet.IUPAC.IUPACAmbiguousDNA())
        newRecord = Bio.SeqRecord.SeqRecord(seq_reformatted,id=record.id, name=record.name,description=primerDescription)
        return newRecord

def plate_vertical_primer_order_fasta(record_iterator,plasmid_name,prefix,startindex,startwell,platetype="numerical"):
    plate = IDT_plate_object(prefix,startindex,startwell,platetype,ignoreWellOverlap=False)
    for record in record_iterator:

        plasmid = plasmid_object(plasmid_name)  ##Needs to know which plasmid is being cloned into    
        ##Write a plasmid file
        plasmid.writePlasmidFile(record,plate=plate) ##Inserts CDS into plasmid.  Uses plate for description type stuff.  Plate optional.

        print("________________________________________________________")
        #if record.seq[:45].translate()[0] != "M" or record.seq[-45:].translate()[14] != "*":
        #    print "***",record.id ,"*** not properly translated"
        #else:
        print("Record id:",record.id)
        print("Record description:",record.description)
        print("")
        if "GGTCTC" in record.seq or "CCAGAG" in record.seq:
            print("Incompatible with BsaI")

            print("GGTCTC:",record.seq.find("GGTCTC"))
            print("CCAGAG:",record.seq.find("CCAGAG"))
        if "CGTCTC" in record.seq or "GCAGAG" in record.seq:
            print("Incompatible with BsmBI")    
            print("CGTCTC:",record.seq.find("CGTCTC"))
            print("GCAGAG:",record.seq.find("GCAGAG"))
        if "GCGGCCGC" in record.seq:
            print("Incompatible with NotI")  
            print("GCGGCCGC:",record.seq.find("GCGGCCGC"))

        ###Related primer pairs in a vertical orientation.  E.g. A1,B1 (numeric incrementing) or A1,A2 (alphabet incremeting)
        #plate.newPrimer(record)
        forwardPrimer = plasmid.formatGibsonPrimerForward(record,24)
        reversePrimer = plasmid.formatGibsonPrimerReverse(record,24)
        plate.newPrimerPairVertical(forwardPrimer,reversePrimer)

        ###Ordering only single primers at a time.
        #plate.newPrimer(record)
        
    plate.df.to_excel("IDT.xls",index=False) 
    
def single_primer_order_fasta(record_iterator,plasmid_name,prefix,startindex,startwell,platetype="numerical"):
    plate = IDT_plate_object(prefix,startindex,startwell,platetype,ignoreWellOverlap=True)
    for record in record_iterator:

        plasmid = plasmid_object(plasmid_name)  ##Needs to know which plasmid is being cloned into    
        ##Write a plasmid file
        plasmid.writePlasmidFile(record,plate=plate) ##Inserts CDS into plasmid.  Uses plate for description type stuff.  Plate optional.

        print("________________________________________________________")
        #if record.seq[:45].translate()[0] != "M" or record.seq[-45:].translate()[14] != "*":
        #    print "***",record.id ,"*** not properly translated"
        print("Record id:",record.id)
        print("Record description:",record.description)
        print("")
        if "GGTCTC" in record.seq or "CCAGAG" in record.seq:
            print("Incompatible with BsaI")

            print("GGTCTC:",record.seq.find("GGTCTC"))
            print("CCAGAG:",record.seq.find("CCAGAG"))
        if "CGTCTC" in record.seq or "GCAGAG" in record.seq:
            print("Incompatible with BsmBI")    
            print("CGTCTC:",record.seq.find("CGTCTC"))
            print("GCAGAG:",record.seq.find("GCAGAG"))
        if "GCGGCCGC" in record.seq:
            print("Incompatible with NotI")    
            print("GCGGCCGC:",record.seq.find("GCGGCCGC"))

        forwardPrimer = plasmid.formatGibsonPrimerForward(record,24)
        reversePrimer = plasmid.formatGibsonPrimerReverse(record,24)
        plate.newPrimer(forwardPrimer)
        plate.newPrimer(reversePrimer)
        
    return plate.getCSV() 
    





if __name__ == "__main__":
    CDS_fasta_file = "photinus_pyralyis_151_v1_Trinity_cleaned.fasta.transdecoder.cds_annotated"
    handle = open(CDS_fasta_file,"rb")
    record_iterator = list(Bio.SeqIO.parse(fasta_handle,"fasta"))
    plate_vertical_primer_order_fasta(record_iterator,"pHis8_4","TRF_x",1,"A1")
else:
    print("importing auto_clone")
    
    



