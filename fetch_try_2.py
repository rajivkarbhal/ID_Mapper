import requests
from bs4 import BeautifulSoup
import re
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.detach(), encoding = 'utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.detach(), encoding = 'utf-8')

# Its demo change file
# Second change

# Another change

class gene:
    def __init__(self,unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl):
        #self.gnurl = gnurl
        self.unurl = unurl
        self.strurl = strurl
        self.ensurl = ensurl
        self.ucscurl = ucscurl
        self.omimurl = omimurl
        self.snpurl = snpurl
        self.varurl = varurl
        self.hlmurl = hlmurl
        self.cddurl = cddurl
        self.bgrdurl = bgrdurl
        self.minturl = minturl
        self.stringurl = stringurl
        self.cprturl = cprturl
        self.inturl = inturl
        self.keggurl = keggurl
        self.gobpurl = gobpurl
        self.goccurl = goccurl
        self.gomfurl = gomfurl

    def gnsmb(self):
        gnid = []
        res = []
        gnsnpn=""
        try:
            urlncbi = self.gnurl
            response = requests.get(urlncbi)
            cot = response.text
            blocks = response.text.split("\n")
            ctst=[]
            for item in enumerate(blocks):
                if (pnm := re.match(r".+\<Count\>(\d+)\<\/Count\>.+", str(item), re.IGNORECASE)):
                    gncount = pnm.group(1)
            #print ("Gene Count: ",  gncount)
            for item in enumerate(blocks): 
                if (pnm := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(item), re.IGNORECASE)):
                    acc = pnm.group(1)
                    gnid.append(acc)
            for ids in (gnid): 
                gnsnpn="" 
                fetchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s"%ids
                NCBI_response = requests.get(fetchurl)
                NCBI_cot = NCBI_response.text
                NCBI_blocks = NCBI_response.text.split("},")                
                gnsnpn = ids 
                for NCBI_item in enumerate(NCBI_blocks):                   
                    if (pnm2 := re.match(r".+taxname \"(.+)\"\,.+", str(NCBI_item), re.IGNORECASE)):
                        org = pnm2.group(1)
                        org = re.sub(r'\".+', '', org)
                        gnsnpn = gnsnpn + "_" + org 
                    if (pnm2 := re.match(r".+locus \"(.\S+)\"\,.+", str(NCBI_item), re.IGNORECASE)):
                        Gnnm = pnm2.group(1)
                        gnsnpn = gnsnpn + "_" + Gnnm
                        res.append(gnsnpn)
            return res
            NCBI_blocks.clear()
            blocks.clear()       
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def gbendd(self):
        idall = []
        try:
            urlged = self.gnurl
            ged_response = requests.get(urlged)
            print (ged_response.text)
            ged_blocks = ged_response.text.split("\n")
            for ged_item in enumerate(ged_blocks):
                if (pnm2 := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(ged_item), re.IGNORECASE)):
                    allid = pnm2.group(1)
                    if allid != selid:
                        idall.append(allid)
            result = ",".join(idall)
            mapacurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=acc&retmode=text"%result
            mapac_response = requests.get(mapacurl)
            mapac = mapac_response.text
            allac = mapac.split("\n")
            return allac
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def unimap (self):
        idall = []
        try:
            urlunmap = self.unurl
            unp_response = requests.get(urlunmap)
            maunipac = unp_response.text
            #print (maunipac)
            alluniac = maunipac.split("\n")
            alluniac.pop(0)
            return alluniac
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def pdbmap (self):
        stridall = []
        try:
            strmap = self.strurl
            #print (strmap)
            str_response = requests.get(strmap)
            mastrac = str_response.text
            pdbcot = mastrac.split("\n")
            one = "https://data.rcsb.org/graphql?query=%7Bentry(entry_id:"
            two = ")%7Bpolymer_entities%7Brcsb_entity_source_organism%7Bncbi_scientific_name%7D%7D%7D%7D"
            for strcd in enumerate(pdbcot):
                 if (regpdb := re.match(r".+\"identifier\" \: \"(\w+)\"\,.+", str(strcd), re.IGNORECASE)):          
                    pdbid = regpdb.group(1)
                    pdbid2 = "\"" + pdbid + "\""
                    comb = one + pdbid2 + two
                    pdbac_response = requests.get(comb)
                    pdcnt= pdbac_response.text
                    if (pdso := re.match(r".+\"ncbi_scientific_name\":\"(.+)\"\}.+", str(pdcnt), re.IGNORECASE)):
                        org = pdso.group(1)
                        comp = pdbid + " - " + org
                        stridall.append(comp)
            return stridall
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def ensmap (self):        
        try:
            ensidall = []
            ensdb = self.ensurl 
            response = requests.get(ensdb)
            blocks = response.text.split("\n")
            for gnidc in enumerate(blocks):
            	if (pnm := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(gnidc), re.IGNORECASE)):
            		bul = pnm.group(1)
            		gnfetchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s"%bul
            		gnfet_response = requests.get(gnfetchurl)
            		gnfet_blocks = gnfet_response.text.split("},")
            		for gnitem in enumerate(gnfet_blocks):
            			if (srcorg := re.match(r".+taxname \"(.+)\"\,", str(gnitem), re.IGNORECASE)):
            				orgnm = srcorg.group(1)
            				src_tx = re.sub("\".+", '', orgnm)
            			if (gnpnm := re.match(r".+db \"Ensembl\"\,", str(gnitem), re.IGNORECASE)):
            				if (gnidm := re.match(r".+tag str \"(\S+)\"", str(gnitem), re.IGNORECASE)):
            					ensmid = gnidm.group(1)
            					comb = ensmid + " - " + src_tx
            					ensidall.append(comb)
            return ensidall
            gnfet_blocks.clear()
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def ucscmap (self):        
        try:
            ucscidall = []
            ucscdb = self.ucscurl          
            response = requests.get(ucscdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (pnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		ucsorg = pnm.group(1)
            		bul = pnm.group(2)
            		ucscomb = bul + " - " + ucsorg
            		ucscidall.append(ucscomb)            
            return ucscidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def omimmap (self):        
        try:
            omimidall = []
            omimdb = self.omimurl 
            response = requests.get(omimdb)
            blocks = response.text.split("\n")
            for gnidc in enumerate(blocks):
            	if (pnm := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(gnidc), re.IGNORECASE)):
            		bul = pnm.group(1)
            		omfetchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=omim&id=%s"%bul
            		omfet_response = requests.get(omfetchurl)
            		omfet_blocks = omfet_response.text.split("\n")
            		for line in omfet_blocks:
            			if (ommt := re.match(r".+<LinkSetDb>", str(line), re.IGNORECASE)):
            				for eln in omfet_blocks:
            					if (ompnm := re.match(r".+\<Id\>(\d+)\<\/Id\>", str(eln), re.IGNORECASE)):
            						ommids = ompnm.group(1)
            						omimidall.append(ommids)
            omimidall.pop(0)					
            return omimidall
            omfet_blocks.clear()
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def snpmap (self):        
        try:
            snpidall = []
            snpdb = self.snpurl
            response = requests.get(snpdb)
            blocks = response.text.split("\n")
            for gnidc in enumerate(blocks):
            	if (pnm := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(gnidc), re.IGNORECASE)):
            		bul = pnm.group(1)
            		snpfetchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=snp&id=%s"%bul
            		snpfet_response = requests.get(snpfetchurl)
            		snpfet_blocks = snpfet_response.text.split("\n")
            		for line in snpfet_blocks:
            			if (snpmt := re.match(r".+<LinkSetDb>", str(line), re.IGNORECASE)):
            				for eln in snpfet_blocks:
            					if (snppnm := re.match(r".+\<Id\>(\d+)\<\/Id\>", str(eln), re.IGNORECASE)):
            						snpids = snppnm.group(1)
            						#print (snpids)
            						snpidall.append(snpids)
            snpidall.pop(0)					
            return snpidall
            snpfet_blocks.clear()
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def varmap (self):        
        try:
            varidall = []
            vardb = self.varurl
            response = requests.get(vardb)
            blocks = response.text.split("\n")
            for gnidc in enumerate(blocks):
            	if (pnm := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(gnidc), re.IGNORECASE)):
            		bul = pnm.group(1)
            		varfetchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=dbvar&id=%s"%bul
            		varfet_response = requests.get(varfetchurl)
            		varfet_blocks = varfet_response.text.split("\n")
            		for line in varfet_blocks:
            			if (varmt := re.match(r".+<LinkSetDb>", str(line), re.IGNORECASE)):
            				for eln in varfet_blocks:
            					if (varpnm := re.match(r".+\<Id\>(\d+)\<\/Id\>", str(eln), re.IGNORECASE)):
            						varids = varpnm.group(1)
            						#print (snpids)
            						varidall.append(varids)
            varidall.pop(0)					
            return varidall
            varfet_blocks.clear()
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def hmlmap (self):        
        try:
            hlmidall = []
            hlmdb = self.hlmurl
            response = requests.get(hlmdb)
            blocks = response.text.split("\n")
            for gnidc in enumerate(blocks):
            	if (pnm := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(gnidc), re.IGNORECASE)):
            		bul = pnm.group(1)
            		hlmfetchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=homologene&id=%s"%bul
            		#print (hlmfetchurl)
            		hlmfet_response = requests.get(hlmfetchurl)
            		#print (hlmfet_response.text)
            		hlmfet_blocks = hlmfet_response.text.split("\n")
            		for line in hlmfet_blocks:
            			if (hlmmt := re.match(r".+<LinkSetDb>", str(line), re.IGNORECASE)):
            				for eln in hlmfet_blocks:
            					if (hlmpnm := re.match(r".+\<Id\>(\d+)\<\/Id\>", str(eln), re.IGNORECASE)):
            						hlmids = hlmpnm.group(1)
            						#print (snpids)
            						hlmidall.append(hlmids)
            #hlmidall.pop(0)					
            return hlmidall
            hlmfet_blocks.clear()
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def cddmap (self):        
        try:
            cddidall = []
            cdddb = self.cddurl
            response = requests.get(cdddb)
            blocks = response.text.split("\n")
            for gnidc in enumerate(blocks):
            	if (pnm := re.match(r".+\<Id\>(\d+)\<\/Id\>.+", str(gnidc), re.IGNORECASE)):
            		bul = pnm.group(1)
            		cddfetchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=cdd&id=%s"%bul
            		#print (hlmfetchurl)
            		cddfet_response = requests.get(cddfetchurl)
            		#print (hlmfet_response.text)
            		cddfet_blocks = cddfet_response.text.split("\n")
            		for line in cddfet_blocks:
            			if (cddmt := re.match(r".+<LinkSetDb>", str(line), re.IGNORECASE)):
            				for eln in cddfet_blocks:
            					if (cddpnm := re.match(r".+\<Id\>(\d+)\<\/Id\>", str(eln), re.IGNORECASE)):
            						cddids = cddpnm.group(1)
            						#print (snpids)
            						cddidall.append(cddids)
            cddidall.pop(0)					
            return cddidall
            cddfet_blocks.clear()
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def biogridmap (self):        
        try:
            bgrdidall = []
            bgrddb = self.bgrdurl          
            response = requests.get(bgrddb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (bgrdpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		bgrdorg = bgrdpnm.group(1)
            		bul = bgrdpnm.group(2)
            		bgrdcomb = bul + " - " + bgrdorg
            		bgrdidall.append(bgrdcomb)            
            return bgrdidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def mintidmap (self):        
        try:
            mintidall = []
            mintdb = self.minturl          
            response = requests.get(mintdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (mintpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		mintorg = mintpnm.group(1)
            		bul = mintpnm.group(2)
            		mintcomb = bul + " - " + mintorg
            		mintidall.append(mintcomb)            
            return mintidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def stringidmap (self):        
        try:
            stringidall = []
            stringdb = self.stringurl
            response = requests.get(stringdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (stringpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		stringorg = stringpnm.group(1)
            		bul = stringpnm.group(2)
            		stringcomb = bul + " - " + stringorg
            		stringidall.append(stringcomb)            
            return stringidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def cmpotmap (self):        
        try:
            cprtidall = []
            cprtdb = self.cprturl
            response = requests.get(cprtdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (cprtpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		cprtorg = cprtpnm.group(1)
            		bul = cprtpnm.group(2)
            		cprtcomb = bul + " - " + cprtorg
            		cprtidall.append(cprtcomb)            
            return cprtidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def intactmap (self):        
        try:
            intactidall = []
            intactdb = self.inturl
            response = requests.get(intactdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (intactpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		intactorg = intactpnm.group(1)
            		bul = intactpnm.group(2)
            		intactcomb = bul + " - " + intactorg
            		intactidall.append(intactcomb)            
            return intactidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def keggmap (self):        
        try:
            keggidall = []
            keggdb = self.keggurl
            response = requests.get(keggdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (keggpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		keggorg = keggpnm.group(1)
            		bul = keggpnm.group(2)
            		keggcomb = bul + " - " + keggorg
            		keggidall.append(keggcomb)            
            return keggidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def gobpmap (self):        
        try:
            gobpidall = []
            gobpdb = self.gobpurl
            response = requests.get(gobpdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (gobppnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		gobporg = gobppnm.group(1)
            		bul = gobppnm.group(2)
            		gobpcomb = gobporg + " - " + bul
            		gobpidall.append(gobpcomb)            
            return gobpidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def goccmap (self):        
        try:
            goccidall = []
            goccdb = self.goccurl
            response = requests.get(goccdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (goccpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		goccorg = goccpnm.group(1)
            		bul = goccpnm.group(2)
            		gocccomb = goccorg + " - " + bul
            		goccidall.append(gocccomb)            
            return goccidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")

    def gomfmap (self):        
        try:
            gomfidall = []
            gomfdb = self.gomfurl
            response = requests.get(gomfdb)
            blocks = response.text.split("\n")
            blocks.pop(0)       
            for gnidc in blocks:
            	if (gomfpnm := re.match(r"\S+\t(.+)\t(.+)\;", str(gnidc), re.IGNORECASE)):
            		gomforg = gomfpnm.group(1)
            		bul = gomfpnm.group(2)
            		gomfcomb = gomforg + " - " + bul
            		gomfidall.append(gomfcomb)            
            return gomfidall
            blocks.clear()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")          




#### *********Take input from GUI********* ###################

#variable_received = sys.argv[1:]
#gen = variable_received[0]
#dblist = variable_received[1:]
availgn = []
gacid = []
unurl=""
strurl=""
ensurl=""
ucscurl=""
omimurl=""
snpurl=""
varurl=""
hlmurl=""
cddurl=""
bgrdurl=""
minturl=""
stringurl=""
cprturl=""
inturl=""
keggurl=""
gobpurl=""
goccurl=""
gomfurl=""
#############################################

gen = sys.argv[1]
input_mode = sys.argv[1:]
dblist = input_mode[1:]
dblt = len(dblist)
if dblt == 1:
	print ("Please select database/s");
else:
	del dblist[-1]
	for ecdb in dblist:
		if ecdb == "GenBanl/ENA/DDBJ":		
			for abc in gacid:
				mapgnacid=[]
				spldt = []
				acno = []
				acunq = []
				#print (abc,"\n")
				spldt = abc.split('/')
				fel = spldt[0]
				if (gmt := re.match(r".+\_(\w+)$", str(fel), re.IGNORECASE)):
					symgen = gmt.group(1)
					print (symgen,"\n")
				for rko in spldt:
					#print (rko,"\n")
					if (mtch := re.match(r"(\d+)\_.+", str(rko), re.IGNORECASE)):
						ecgnid = mtch.group(1)
						acurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=nuccore&id=%s"%ecgnid
						#print (acurl)
						geneacmap = gene(acurl)
						acno=geneacmap.gbendd()
				spldt.clear()
				for uac in acno:
					mapgnacid.append(uac)
				acunq = list(set(mapgnacid))
				for tej in acunq:
					print (tej)
				acno.clear()
				mapgnacid.clear()
				acunq.clear()

	#Calling function for UniProt mapping accession	
		if ecdb == "UniProtKB/Swiss-Prot":
			print ("********** Mapped ID on UniProtKB/Swiss-Prot database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					prturl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Corganism_name&format=tsv&query=(gene%3A{})+AND+(reviewed:true)".format(ecg)
					geneid = gene(prturl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					prtacno=geneid.unimap()
					length = len(prtacno)
					if length == 1:
						print(ecg, "is not mapped on UniProtKB/Swiss-Prot database")
					else:
						for unitem in (prtacno):
							if (pnm2 := re.match(r"(\w+)\t(\w+)\t(.+)", str(unitem), re.IGNORECASE)):
								unpac = pnm2.group(1)
								unpacorg = pnm2.group(3)
								uncomb = unpac + " - " + unpacorg
								print (uncomb) 
			else:
				print ("\nGene: ", gen)
				prturl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Corganism_name&format=tsv&query=(gene%3A{})+AND+(reviewed:true)".format(gen)
				geneid = gene(prturl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				prtacno=geneid.unimap()
				length = len(prtacno)
				if length == 1:
					print(gen, "is not mapped on UniProtKB/Swiss-Prot database")
				else:
					for unitem in (prtacno):
						if (pnm2 := re.match(r"(\w+)\t(\w+)\t(.+)", str(unitem), re.IGNORECASE)):
							unpac = pnm2.group(1)
							unpacorg = pnm2.group(3)
							uncomb = unpac + " - " + unpacorg
							print (uncomb)

	#Calling function for PDB mapping	
		if ecdb == "PDB":
			print ("\n********** Mapped ID on PDB database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					pdburl = "https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%22query%22:%7B%22type%22:%22terminal%22,%22label%22:%22text%22,%22service%22:%22text%22,%22parameters%22:%7B%22attribute%22:%22rcsb_entity_source_organism.rcsb_gene_name.value%22,%22operator%22:%22exact_match%22,%22negation%22:false,%22value%22:%22{}%22%7D%7D,%22return_type%22:%22entry%22,%22request_options%22:%7B%22paginate%22:%7B%22start%22:0,%22rows%22:500%7D,%22results_content_type%22:%5B%22experimental%22%5D,%22sort%22:%5B%7B%22sort_by%22:%22rcsb_entry_info.resolution_combined%22,%22direction%22:%22asc%22%7D%5D,%22scoring_strategy%22:%22combined%22%7D%7D".format(ecg)
					geneid = gene(unurl,pdburl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					stracno=geneid.pdbmap()
					length = len(stracno)
					if length == 0:
						print(ecg, "is not mapped on PDB database")
					for stritem in (stracno):
						print (stritem)
			else:
				print ("\nGene: ", gen)
				pdburl = "https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%22query%22:%7B%22type%22:%22terminal%22,%22label%22:%22text%22,%22service%22:%22text%22,%22parameters%22:%7B%22attribute%22:%22rcsb_entity_source_organism.rcsb_gene_name.value%22,%22operator%22:%22exact_match%22,%22negation%22:false,%22value%22:%22{}%22%7D%7D,%22return_type%22:%22entry%22,%22request_options%22:%7B%22paginate%22:%7B%22start%22:0,%22rows%22:500%7D,%22results_content_type%22:%5B%22experimental%22%5D,%22sort%22:%5B%7B%22sort_by%22:%22rcsb_entry_info.resolution_combined%22,%22direction%22:%22asc%22%7D%5D,%22scoring_strategy%22:%22combined%22%7D%7D".format(gen)
				geneid = gene(unurl,pdburl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				stracno=geneid.pdbmap()
				length = len(stracno)
				if length == 0:
					print(gen, "is not mapped on PDB database")
				else:
					for stritem in (stracno):
						print (stritem)

	#Calling function for ENSEMBL mapping	
		if ecdb == "ENSEMBL":
			print ("\n********** Mapped ID on ENSEMBL database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					gnnmurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % ecg
					geneid = gene(unurl,strurl,gnnmurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.ensmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on ENSEMBL database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				gnnmurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % gen
				geneid = gene(unurl,strurl,gnnmurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.ensmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on ENSEMBL database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for UCSC mapping	
		if ecdb == "UCSC":
			print ("\n********** Mapped ID on UCSC database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					ucurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_ucsc&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.ucscmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on UCSC database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				ucurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_ucsc&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				geneid = gene(unurl,strurl,ensurl,ucurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.ucscmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on UCSC database")
				else:
					for ensitem in (ensacno):
						print (ensitem)					

	#Calling function for OMIM mapping	
		if ecdb == "OMIM":
			print ("\n********** Mapped ID on OMIM database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					omurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % ecg
					geneid = gene(unurl,strurl,ensurl,ucscurl,omurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno = geneid.omimmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on OMIM database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				omurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % gen
				geneid = gene(unurl,strurl,ensurl,ucscurl,omurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.omimmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on OMIM database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for SNP mapping	
		if ecdb == "SNP":
			print ("\n********** Mapped ID on SNP database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					spurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % ecg
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,spurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno = geneid.snpmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on SNP database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				spurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % gen
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,spurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.snpmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on SNP database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for dbVar mapping	
		if ecdb == "dbVar":
			print ("\n********** Mapped ID on dbVar database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					vrurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % ecg
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,vrurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno = geneid.varmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on dbVar database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				vrurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % gen
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,vrurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.varmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on dbVar database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for dbVar mapping	
		if ecdb == "HomoloGene":
			print ("\n********** Mapped ID on HomoloGene database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					hlurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % ecg
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno = geneid.hmlmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on HomoloGene database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				hlurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % gen
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.hmlmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on HomoloGene database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for CDD mapping	
		if ecdb == "CDD":
			print ("\n********** Mapped ID on CDD database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					cdurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % ecg
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cdurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno = geneid.cddmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on CDD database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				cdurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]&retmax=100000" % gen
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cdurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.cddmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on CDD database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for BioGRID mapping	
		if ecdb == "BioGRID":
			print ("\n********** Mapped ID on BioGRID database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					bdurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_biogrid&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.biogridmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on BioGRID database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				bdurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_biogrid&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.biogridmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on BioGRID database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for MINT mapping	
		if ecdb == "MINT":
			print ("\n********** Mapped ID on MINT database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					mturl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_mint&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,mturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.mintidmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on MINT database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				mturl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_mint&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,mturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.mintidmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on MINT database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for STRING mapping	
		if ecdb == "STRING":
			print ("\n********** Mapped ID on STRING database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					sturl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_string&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,sturl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.stringidmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on STRING database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				sturl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_string&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				#print (sturl)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,sturl,cprturl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.stringidmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on STRING database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for Complex Portal mapping	
		if ecdb == "Complex Portal":
			print ("\n********** Mapped ID on Complex Portal database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					cpurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_complexportal&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cpurl,inturl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.cmpotmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on Complex Portal database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				cpurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_complexportal&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				#print (sturl)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cpurl,inturl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.cmpotmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on Complex Portal database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for IntAct database mapping	
		if ecdb == "IntAct":
			print ("\n********** Mapped ID on IntAct database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					iaurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_intact&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,iaurl,keggurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.intactmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on Complex Portal database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				iaurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_intact&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				#print (sturl)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,iaurl,keggurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.intactmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on Complex Portal database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for KEGG database mapping	
		if ecdb == "KEGG":
			print ("\n********** Mapped ID on KEGG database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					kgurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_kegg&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,kgurl,gobpurl,goccurl,gomfurl)
					ensacno=geneid.keggmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on Complex Portal database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				kgurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cxref_kegg&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				#print (sturl)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,kgurl,gobpurl,goccurl,gomfurl)
				ensacno=geneid.keggmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on Complex Portal database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for GOBP database mapping	
		if ecdb == "Gene Ontology (Biological_process)":
			print ("\n********** Mapped ID on Gene Ontology (Biological_process) database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					gbpurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cgo_p&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gbpurl,goccurl,gomfurl)
					ensacno=geneid.gobpmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on Gene Ontology (Biological_process) database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				gbpurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cgo_p&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				#print (sturl)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gbpurl,goccurl,gomfurl)
				ensacno=geneid.gobpmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on Gene Ontology (Biological_process) database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for GOBP database mapping	
		if ecdb == "Gene Ontology (Cellular_component)":
			print ("\n********** Mapped ID on Gene Ontology (Cellular_component) database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					gcurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cgo_c&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,gcurl,gomfurl)
					ensacno=geneid.goccmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on Gene Ontology (Cellular_component) database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				gcurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cgo_c&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				#print (sturl)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,gcurl,gomfurl)
				ensacno=geneid.goccmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on Gene Ontology (Cellular_component) database")
				else:
					for ensitem in (ensacno):
						print (ensitem)

	#Calling function for GOBP database mapping	
		if ecdb == "Gene Ontology (Molecular_function)":
			print ("\n********** Mapped ID on Gene Ontology (Molecular_function) database **********")
			if (sep := re.match(r".+\,.+", str(gen), re.IGNORECASE)):
				gnl = gen.split(",")
				for ecg in gnl:
					print ("\nGene: ", ecg)
					gmurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cgo_f&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(ecg)
					geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gmurl)
					ensacno=geneid.gomfmap()
					length = len(ensacno)
					if length == 0:
						print(ecg, "is not mapped on Gene Ontology (Molecular_function) database")
					else:
						for ensitem in (ensacno):
							print (ensitem) 
			else:
				print ("\nGene: ", gen)
				gmurl = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Corganism_name%2Cgo_f&format=tsv&query=%28%28gene%3A{}%29%29+AND+%28reviewed%3Atrue%29".format(gen)
				#print (sturl)
				geneid = gene(unurl,strurl,ensurl,ucscurl,omimurl,snpurl,varurl,hlmurl,cddurl,bgrdurl,minturl,stringurl,cprturl,inturl,keggurl,gobpurl,goccurl,gmurl)
				ensacno=geneid.gomfmap()
				length = len(ensacno)
				if length == 0:
					print(gen, "is not mapped on Gene Ontology (Molecular_function) database")
				else:
					for ensitem in (ensacno):
						print (ensitem)																																												
