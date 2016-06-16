#! /usr/bin/python
##Timothy R. Fallon 2016
##Gibson cloning primer producer.


import tornado
import tornado.ioloop
import tornado.web
import os, uuid
import os.path
import StringIO
import Bio
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.Alphabet
import glob
import re
import logging
import base64
import auto_clone
import zipfile
 
TCP_PORT = 5686
 
class Userform(tornado.web.RequestHandler):
    def get(self):
        self.render("fileuploadform.html")
 
 
class Upload(tornado.web.RequestHandler):
    def post(self):
    
        os.chdir("/lab/solexa_weng/Seq_data/Tim/Scripts/auto_clone")
    
        self.write("<html><font face=\"courier\">")
        self.write("<div style =\"width: 5000px; overflow:scroll\">") ##Some tricks to have scrolling text.
        plasmid_name = self.get_argument('plasmid_name', '')
        fasta_text = self.get_argument('fasta_text','').strip().lower().replace('\n','').replace(' ','')
        record_id = self.get_argument('record_id','').strip()
        record_name = self.get_argument('record_name','').strip()
        record_description = self.get_argument('record_description','').strip()
        
        primer_prefix = self.get_argument('primer_prefix','').strip()
        primer_index = int(self.get_argument('primer_index','').strip())
	starting_well = self.get_argument('starting_well','').strip()
        
        IDT_format = self.get_argument('IDT_format','').strip().strip('\n')
        
        ##Randomly assign name for temporary directory
        tmp_str = base64.b64encode(os.urandom(16),['-','_'])
        tmp_dir_name = "./output_folders/"+tmp_str
        os.mkdir(tmp_dir_name)
        current_dir = os.getcwd()
    
        ##Setup input stream.
        if fasta_text == '':
            ##Text box is empty, is there a file?
            try:
                fileinfo = self.request.files['filearg'][0]
                fname = fileinfo['filename']
                extn = os.path.splitext(fname)[1]
                file_data = fileinfo['body'].strip()
		file_data = str(file_data.replace("\r","\n"))
		input_stream = StringIO.StringIO(file_data)
            except:
                self.finish("No file uploaded or sequence input.")
                return
        else:
            fasta_text = fasta_text.strip()
            fasta_text = str(fasta_text.replace("\r","\n"))
            input_stream = StringIO.StringIO(fasta_text)  
    
        ##Parse input stream. 
        if fasta_text == '':
            ##A file was uploaded. Try to parse it as fasta.
            record_iterator = list(Bio.SeqIO.parse(input_stream,"fasta"))
        else:
            if fasta_text[0] == ">":
                ##Seems like the text box has fasta, parse as such.
                record_iterator = list(Bio.SeqIO.parse(input_stream,"fasta"))
		print input_stream.getvalue()
            elif fasta_text[0] == 'a' or fasta_text[0] == 'c' or fasta_text[0] == 't' or fasta_text[0] == 'g':
                ##Seems like raw sequence, convert to record.
                record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(fasta_text\
                ,Bio.Alphabet.IUPAC.IUPACAmbiguousDNA()),id=record_id, name=record_name
                ,description=record_description)
                record_iterator = [record]
            else:
                self.finish("couldn't parse FASTA text box")
                return
    
	self.write("Loaded "+str(len(record_iterator))+" record(s) successfully")        
        self.write("<br>Ordering gibson primers for plasmid: "+plasmid_name)
        self.write("<br>IDT format is: "+IDT_format)

        os.chdir(tmp_dir_name)
        if IDT_format == "plate":
            auto_clone.plate_vertical_primer_order_fasta(record_iterator,plasmid_name,primer_prefix,primer_index,starting_well)
            self.write("<br>Produced IDT primer order file:")
            self.write("<br><a href="+"output_folders/"+tmp_str+"/IDT.xls>IDT.xls</a><br>")
	    self.write("<br>Order from IDT here <a href=\"https://www.idtdna.com/site/order/plate/index/dna/1800\">https://www.idtdna.com/site/order/plate/index/dna/1800</a>")            
        elif IDT_format == "single":
	    self.write("<br>Order from IDT here <a href=\"https://www.idtdna.com/site/order/oligoentry\">https://www.idtdna.com/site/order/oligoentry</a>")            
            self.write("<br>Copy paste below into IDT (CSV bulk input):<br><br>")
            result = auto_clone.single_primer_order_fasta(record_iterator,plasmid_name,primer_prefix,primer_index,starting_well) 
            print result
            result = result.replace('\n','<br>')
            self.write(result)
        os.chdir(current_dir)
        
        self.write("<br>Produced plasmid files:")
        zipname = tmp_dir_name+"/plasmids.zip"
        zip_handle = None
        for root, subFolders, files in os.walk(tmp_dir_name):
            for filename in files:
                if filename.lower().endswith(".gb"):
                    self.write("<br>"+filename)
                    if zip_handle == None: 
                        zip_handle = zipfile.ZipFile(zipname,"w") ##Open the zip file.
                    zip_handle.write(tmp_dir_name+"/"+filename,filename) ##Append files to the zip file.
        
        if not zip_handle == None:
            zip_handle.close()
        
        if os.path.isfile(zipname):
        ##Check if the file was actually produced
            self.write("<br><a href="+"output_folders/"+tmp_str+"/plasmids.zip>Plasmid .zip download</a>")
        else:
            self.write("<br>Couldn't produce plasmid files (likely due to missing reference .gb internally)")
    
    
        self.write("</div>")

        ##Making the return text.  Tries to return only unique records (non-matching SeqID & sequence)
        ##Since it holds the sequences in memory to compare if duplicates exist, this could deal badly with large memory fasta files.
        self.finish("</font></html>")
        #self.finish(cname + " is uploaded!! Check %s folder" %__UPLOADS__)
        return 0
 
 
application = tornado.web.Application([
        (r"/", Userform),
        (r"/upload", Upload),
        (r"/output_folders/(.*)",tornado.web.StaticFileHandler,{'path':"./output_folders/"}) ##http://stackoverflow.com/questions/10165665/using-tornado-how-do-i-serve-static-files-and-serve-a-favicon-ico-from-a-differ
        ], debug=True)
 
 
if __name__ == "__main__":
    
    pid = str(os.getpid())
    f = open('fasta_lookup_server.pid', 'w')
    f.write(pid)
    f.close()

    application.listen(TCP_PORT)
    tornado.ioloop.IOLoop.instance().start()
