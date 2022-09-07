#Importing BeautifulSoup
from bs4 import BeautifulSoup
import requests as rq
import batch_splitter
#Variable used to assign a UID to each study
uid=0

if __name__ == "__main__":   
    #Generating the search URL
    base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    search_url = base + 'esearch.fcgi?db=popset&term=aves+AND+Population+Study[Filter]&usehistory=y'
    #Sending a get request with the search URL
    search_result = rq.get(search_url)
    
    
    
    #Parsing the result using BeautifulSoup
    doc = BeautifulSoup(search_result.text, "html.parser")
    #Extracting the Query Key
    query_key = doc.querykey.string
    #Extracting the Web Environment
    webenv = doc.webenv.string
    #Extracting total results
    total_results = doc.count.string
    
    
    
    #Fetching the data in batches of 250
    #Calculating needed loops
    loops = int(int(total_results) / 250) + 1
    #Variable that affects where the start point of the group is
    retstart = 0
    #Variable that affects the size of the group
    retmax = 250 
    
    #Creating a loop to extact the data in groups of 250
    for i in range(loops):   
        #Ensures the retmax doesn't exceed the total results
        if (retstart + 250) > int(total_results):
            retmax = int(total_results)
        
        #Generates the fetch URL for the batch
        fetch_url = base + 'efetch.fcgi?db=popset&rettype=fasta&retmode=text'\
            f'&query_key={query_key}&WebEnv={webenv}&retstart={retstart}&retmax={retmax}'
        #Increments the start point to get ready for a new batch
        retstart += 250
        
    
        #Sends a get request for the batch
        fetch_result = rq.get(fetch_url)
        fetch_text = fetch_result.text
        #Uses the study_export function to save each study individually
        batch_splitter.study_export(fetch_text)
        #Increments the UID for the next batch
        uid += 250