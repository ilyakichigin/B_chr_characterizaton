#!/usr/bin/env python

import os
import sys
import pysam
import shlex
import shutil
import subprocess

def run_command(command, 
                verbose=True, dry_run=False, 
                outfile=None, errfile=None,
                stdin=None,
                return_out=False):
    if verbose:
        sys.stderr.write(command)
        if outfile:
            sys.stderr.write(' > %s' % outfile)
        if errfile:
            sys.stderr.write(' 2>> %s' % errfile)
        sys.stderr.write('\n')
    if not dry_run:
        for f in (outfile, errfile):
            if f:
                d = os.path.dirname(f)
                if not os.path.isdir(d):
                    os.makedirs(d)
        command = shlex.split(command)
        p = subprocess.Popen(command,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate(input=stdin)
        if p.returncode != 0:
            raise RuntimeError("%s error:\n%s" % (command[0], err))
        if outfile:
            with open(outfile, 'wb') as o:
                o.write(out)
        if errfile:
            with open(errfile, 'ab') as e:
                e.write(err)
        if return_out:
            return out

def bam_sort(in_file, out_file, name_sort=False, verbose=True, dry_run=False):
    
    if verbose:
        sys.stderr.write('Sorting %s ' % (in_file))
    if name_sort:
        if verbose:
            sys.stderr.write('by read name. Output: %s.\n' % (out_file))
        if not dry_run:
            pysam.sort('-n', '-T', '/tmp/bam_nsort', '-o', out_file, in_file)
    else:
        if verbose:
            sys.stderr.write('by reference coordinate. Output: %s.\n' 
                             % (out_file))
        if not dry_run:
            pysam.sort('-T', '/tmp/bam_sort', '-o', out_file, in_file)

def check_exec(p):
    """Check if executable can be accessed, 
    otherwise raise OSError
    """
    if p.startswith('~'):
        raise OSError(p + ' executable has "~" in path.')
    try:
        fnull = open(os.devnull, 'w')
        subprocess.call(p, stdout=fnull, stderr=fnull)
    except OSError as e:
        raise OSError(p + ' executable not found.')


def check_file(f):
    """Check if executable can be accessed, 
    otherwise raise OSError
    """
    if not os.path.isfile(f):
        raise OSError(f + ' not found.')
    if f.startswith('~'):
        raise OSError(f + ' has "~" in path.')

def copy_to_cwd(source_dir, source_name, dest_name):
    """Copy source file to current folder and name it dest_name
    raise OSError if file exists
    """
    source_path = os.path.join(source_dir, source_name)
    dest_path = os.path.join(os.getcwd(), dest_name)
    if os.path.isfile(dest_path):
        raise OSError(dest_name + ' file exists in current directory.')
    shutil.copy2(source_path, dest_path)