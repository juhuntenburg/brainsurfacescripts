# function to write time and message to log file
def log(log_file, message, logtime=True):
    import time
    with open(log_file, 'a') as f:
        if logtime:
            f.write(time.ctime()+'\n')
        f.write(message+'\n')
        
# function returning generator for sub,hemi tuples fro joblib
def tupler(subjects, hemis):
    for s in subjects:
        for h in hemis:
            yield (s, h)