# matz_tutorial


## Pro Tips:

### TACC

- Always check what directory you're working in before you start using (pwd)
- Every project you start should get it's own directory in SCRATCH, make this using (mkdir)
- Always check the output of a job on TACC by "listing" the contents of your directory. The best way to do this is by (ll -tr), which lists all of your files sorted by most recent. This way you can see what was just created.
- Don't batch more than one job *with the same name* in a row. Check what is currently running using (squeue -u YOURUSERNAME). To cancel a job, if you make an error, find the jobID using squeue -u, then you can (scancel JOBNUMBER)

### File Naming: Underscores are your best friend

- NEVER!! Use spaces in the column names for your data files, instead, use underscores or dashes (e.g., name something "Sampling_Location" instead of "Sampling Location")
- This also goes for filenames, *do not* name your folders, directories, or files with spaces in them. When you read files into software later, this will cause big issues and be a massive pain in your behind.
- *Try* not to use '.' in folder/file names either. This one is less of a rule, and I realize most files end with a ".txt" or ".fasta" extension. However, the '.' character is also a "wildcard" when trying to search for strings. This means that if you use the "grep" command or any other "regrex" command, the '.' character is interpreted to mean ANY character. Again, it's always best to use "_" when trying to separate file names. (e.g., use "My_Sample_ID.fasta" instead of "My.Sample.ID.fasta".
- To make your life easier, come up with a pattern that you use to name files. For example, when I am naming my coral samples, I use a consistent pattern that I will be able to separate later. 
Ex: POR_L21_ADULT.fasta is separated by two "_". Everything before the first "_" gives me the species (this sample is a Porites). Then, I know between the first and second "_", I have the sample id (this is individual L21). Finally, I know it's an adult sample because after the second "_", it says Adult. By naming things with a common pattern, I will be able to "split them on the delimiter" later. A delimiter is just a common character that separates fields - here, I have chosen the underscore ("_"). 
