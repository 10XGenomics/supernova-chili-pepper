import shutil
import subprocess
import tenkit.supernova as tk_sn

def split(args):
    mem_gb = 600
    chunk_defs = [{'__mem_gb': mem_gb}]
    return {'chunks': chunk_defs}

def join(args, outs, chunk_defs, chunk_outs):
    shutil.move(chunk_outs[0].default, outs.default)

def check_exclude(path, ext):
    head=path[:-len(ext)]
    tail=path[-len(ext):]
    if tail != ext:
        raise Exception("file has incorrect extension: " + path)
    return head

def probe_reads(path):
    out=subprocess.check_output(['FastFastbCount','QUIET=True','FASTB='+path]).strip()
    return int(out)

def main(args, outs):
    h1 = check_exclude(args.reads, '.fastb')
    h2 = check_exclude(args.quals, '.qualp')
    h3 = check_exclude(args.bci, '.bci')
    if h1 != h2 or h2 != h3:
        raise Exception( "something wrong with filenames passed in" )

    select_frac=1.0
    if args.downsample is not None:
        if args.downsample.get("target_reads", None) is not None:
            target_nreads = args.downsample["target_reads"]
            actual_nreads = probe_reads( args.reads )

            print "target_nreads={}, actual_nreads={}".format( \
                target_nreads, actual_nreads )

            if target_nreads > actual_nreads:
                select_frac=1.
            else:
                select_frac=1.*target_nreads/actual_nreads

    df_command = ['DF', 'LR_SELECT_FRAC={:f}'.format(select_frac), 'LR='+args.reads,
               'OUT_DIR='+outs.default ]

    if args.pipeline_id is not None:
        df_command.append( "PIPELINE={}".format(args.pipeline_id) )

    if args.known_sample_id is not None:
        df_command.append( "SAMPLE={}".format(args.known_sample_id) )

    if args.addin is not None and "DF" in args.addin:
        df_command.append(args.addin["DF"])

    print " ".join(df_command)

    alarm_bell = tk_sn.SupernovaAlarms(base_dir=outs.default)
    try:
        subprocess.check_call(df_command)
    except subprocess.CalledProcessError:
        ## there was an error of some kind
        ## if any Martian::exit was issued in the C++ code
        ## then we call martian.exit from here
        alarm_bell.post()
        ## if we actually reach this point
        ## it means that Martian::exit() was not called from C++.
        ## issue a martian exit with an unexpected error
        alarm_bell.exit() 
    ## if DF completed successfully,
    ## post any warnings from this stage here
    alarm_bell.post()

