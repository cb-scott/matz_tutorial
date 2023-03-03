# matz_tutorial


## Pro Tips:
- Always check what directory you're working in before you start using (pwd)
- Every project you start should get it's own directory in SCRATCH, make this using (mkdir)
- Always check the output of a job on TACC by "listing" the contents of your directory. The best way to do this is by (ll -tr), which lists all of your files sorted by most recent. This way you can see what was just created.
- Don't batch more than one job *with the same name* in a row. Check what is currently running using (squeue -u YOURUSERNAME). To cancel a job, if you make an error, find the jobID using squeue -u, then you can (scancel JOBNUMBER)
