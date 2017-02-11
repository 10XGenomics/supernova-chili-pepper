import shutil
import subprocess
import os
import tenkit.supernova as tk_sn

def split(args):
    mem_gb = 600
    chunk_defs = [{'__mem_gb': mem_gb}]
    return {'chunks': chunk_defs}

def join(args, outs, chunk_defs, chunk_outs):
    shutil.move(chunk_outs[0].default, outs.default)
    outs.summary_cs = os.path.join(outs.default, "stats", "summary_cs.csv")

def main(args, outs):

    input_dir = os.path.join(args.parent_dir,"a.base")

    cp_command = ['CP','DIR='+input_dir]

    if args.known_sample_id is not None:
        cp_command.append( "SAMPLE={}".format(args.known_sample_id) )
    if args.sample_id is not None:
        cp_command.append( "CS_SAMPLE_ID={}".format(args.sample_id) )
    if args.sample_desc is not None:
        cp_command.append( "CS_SAMPLE_DESC={}".format(args.sample_desc) )

    alarm_bell = tk_sn.SupernovaAlarms(base_dir=args.parent_dir)
    try:
        subprocess.check_call( cp_command )
    except subprocess.CalledProcessError:
        alarm_bell.post()
        ## exit if we haven't already
        alarm_bell.exit()    
    ## post any warnings from this stage here
    alarm_bell.post()

    shutil.move( args.parent_dir, outs.default )
