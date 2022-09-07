#Function used to export each individual study from a batch into a .txt file
def study_export(batch):
    #UID variable used to differentiate between studies
    import aves_sequence_extractor
    #Breaks the results into individual studies
    sequences_list = batch.split('\n\n')
    #Deleting the empty element at the end
    del sequences_list[-1]
    

    #Iterates through each study
    for i in range(len(sequences_list)):
            #Stores the current study data in a new variable
            current_sequence = sequences_list[i]
            
            #Extracting a suitable file name
            file_name = current_sequence.split(' ')
            file_name = '_'.join(file_name[1:4]).replace('/','_')
        
            #Exporting each study into its own text file
            try:
                outfile = open(f"Sequences\{file_name}_UID_{aves_sequence_extractor.uid}.txt", 'w')
                outfile.write(current_sequence)
                outfile.close()
                #Generates the next uid
                aves_sequence_extractor.uid += 1
            #Error returned if a file is skipped
            except:
                print(f"Error on study {file_name}")
    
    print("Batch Complete")
