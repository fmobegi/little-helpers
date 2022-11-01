#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Wrapper for orthoMCL + MySQL management. 
Adapted from Victor de Jager original orthomcl wrapper 2011"""

import sys
import glob
import os
import subprocess
import shlex
import re
import string
import getopt
import MySQLdb

############
# SETTINGS #
############
# template config file
conf_file = '~/orthomclSoftware-v2.0.9/doc/OrthoMCLEngine/Main/orthomcl.config.template'

""" Your SQL login details here!!"""
sql_user = 'user'                         # mysql  root user
sql_passwd = 'pw'     			      # mysql  root pass
bin_dir = '~/orthomclSoftware-v2.0.9/bin/'  # Set your orthomcl bin dir
cpu = 6  # minimum threads to use
chunks = cpu


def rename_IDs(indir, outdir):
    filepart = os.path.split(infile)[1]
    outname = outdir+'/'+filepart


def load_txt(file):
    input = open(file, 'r')
    text = input.read()
    input.close()
    return text


def get_options(options, arguments, longarguments, compulsory):
    """Generic commandline option parsing
    arguments= 'i:'
    longarguments=['cpu= '] # array
    compulsory=['i']
    options={}
    """
    argv = sys.argv[1:]
    print argv
    try:
        opts, leftover = getopt.getopt(argv, arguments, longarguments)
    except getopt.GetoptError as e:
        print e
        usage()
        sys.exit(2)

    if leftover! = []:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        options[opt] = arg

    for element in compulsory:
        if not '-' + element in options.keys():
            usage()
            sys.exit(2)

    return options


def usage():
    sys.stderr.write(
        'Use: %s -i <indir> -n <name_for_db> --cpu=<nr_of_cpus> --chunks=<split fasta file in n chunks, same as cpu if not specified>\n' % sys.argv[0])


def openDB(user, passwd):
    db = MySQLdb.connect(host="localhost", user=user, passwd=passwd)
    cursorDB = db.cursor()
    return cursorDB


def run_command(command):
    command = bin_dir+command
    print "Running: "+command
    result = os.system(command)
    if not result == 0:
        sys.stderr.write("Error running %s\n" % command)
        sys.exit(1)
    return


def run_commands_parallel(commands, cpu, chunks):
    while (commands):
        print "%s commands left to run; " % (len(commands))
        ps = {}
        for process in xrange(cpu):
            try:
                c = commands.pop()
                print "Running: ", c
                args = shlex.split(c)
                # shell should be false, otherwise arguments are not parsed correctly
                subproc = subprocess.Popen(args, shell=False)
                ps[subproc.pid] = subproc
                print "\tWaiting for %d     processes..." % len(ps)
            except:
                pass

        while ps:
            pid, status = os.wait()
            if pid in ps:
                del ps[pid]
            print "\tWaiting for %d processes..." % len(ps)


###############################################################################################################################################################
# main program
arguments = 'n:i:d'
longarguments = ['cpu= ', 'chunks= ']
compulsory = ['n', 'i']
options = {}
## MORE SETTINGS AT THE TOP!! ##

if __name__ == '__main__':
    # parse options
    options = get_options(options, arguments, longarguments, compulsory)
    indir = options['-i']
    name = options['-n']

    # LOCKFILE
    lockname = '/tmp/orthomcl_%s.lock' % name
    print 'Making lockfile'
    if os.path.isfile(lockname):
        sys.stderr.write(
            'Lockfile %s exists, it seems another instance of orthomcl is running on this DB. Please pick another name or remove the lockfile.\n' % lockname)
        sys.exit(1)
    os.system('touch %s' % lockname)

    # SQL PREP
    cursor = openDB(sql_user, sql_passwd)
    if options.has_key('-d'):
        raw_input('About to delete all data from %s (WILL DELETE ALL DATA, CAUTION!)\nPress ctrl-c to abort, any key to continue' % name)
        command = 'DROP SCHEMA IF EXISTS %s;' % name
        print command
        cursor.execute(command)
        try:
            command = 'DROP USER %s;' % name
            print command
            cursor.execute(command)
        except:
            pass
    try:
        command = 'CREATE SCHEMA %s;' % name
        print command
        cursor.execute(command)
    except MySQLdb.ProgrammingError:
        sys.stderr.write(
            'Database "%s" exists, please run with -d to delete and recreate\n' % name)
        print 'Removing lockfile'
        os.system('rm %s' % lockname)
        sys.exit()

    command = "CREATE USER %s IDENTIFIED BY '%s';" % (name, name)
    print command
    cursor.execute(command)
    command = "GRANT ALL PRIVILEGES ON %s.* TO '%s'" % (name, name)
    print command
    cursor.execute(command)
    cursor.close()
    # GET CONF
    text = load_txt(conf_file)
    text = text.replace('dbLogin= ', 'dbLogin=%s' % name)
    text = text.replace('dbPassword= ', "dbPassword=%s" % name)
    text = text.replace('dbConnectString=dbi:mysql:orthomcl:3307',
                        'dbConnectString = dbi:mysql:%s' % name)
    output = open('orthomcl.config', 'w')
    output.write(text)
    output.close()

    # INSTALL SCHEMA
    command = 'orthomclInstallSchema orthomcl.config sql.log'
    run_command(command)

    # CHECK FASTA
    fasta_files = glob.glob(indir+'/*.faa') + \
        glob.glob(indir+'/*.fasta')+glob.glob(indir+'/*.fa')
    print 'Found %i fasta files' % len(fasta_files)
    translate_file = 'IDs_to_tags.txt'
    output_tags = open(translate_file, 'w')
    if not os.path.isdir('adjusted_fasta'):
        os.mkdir('adjusted_fasta')
    for count, filename in enumerate(fasta_files):
        short_filename = os.path.split(filename)[1]
        name_id = str(count+1).zfill(4)
        outname = 'tmp.fasta'
        output = open(outname, 'w')
        orf_count = 1
        text = load_txt(filename).strip().split('\n')
        for line in text:
            if line[0] == '>':
                orf_id = str(orf_count).zfill(4)
                output.write('>%s\n' % orf_id)
                output_tags.write(string.join(
                    [name_id, orf_id, short_filename, line[1:]], '\t')+'\n')
                orf_count += 1
            else:
                output.write(line+'\n')
        output.close()
        command = 'orthomclAdjustFasta %s %s 1' % (name_id, 'tmp.fasta')
        run_command(command)
        os.system('mv %s.fasta adjusted_fasta/%s.fasta' % (name_id, name_id))
    output_tags.close()
    command = 'orthomclFilterFasta adjusted_fasta 20 5'
    run_command(command)

    # RUN BLAST
    if (options.has_key('--cpu')):
        cpu = int(options['--cpu'])

    if (options.has_key('--chunks')):
        chunks = int(options['--chunks'])

    #command = 'formatdb -p T -i goodProteins.fasta'
    command = '/usr/bin/makeblastdb -in goodProteins.fasta -dbtype prot -title goodProteins.fasta -out goodProteins.fasta'
    print command
    os.system(command)

    # split the input file in several files for multi threading
    command = '/usr/bin/fastasplit -f goodProteins.fasta -o ./ -c %s' % chunks
    print command
    os.system(command)

    # get the chunks, but not old .out files
    fastachunks = glob.glob('*chunk*')
    filter = re.compile(".+\d+$", re.IGNORECASE)
    fastachunks = [f for f in fastachunks if filter.search(f)]

    # create the command set
    commandtemplate = '/usr/bin/blastp -db goodProteins.fasta -query %s -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out %s.out'
    commands = [commandtemplate % (chunk, chunk) for chunk in fastachunks]
    outputs = ["%s.out" % chunk for chunk in fastachunks]

    if os.path.isfile('blast.out'):
        print "blastfile already exists"
    else:
        run_commands_parallel(commands, cpu, chunks)
        combinecommand = "cat %s > blast.out" % (' '.join(outputs))
        os.system(combinecommand)
        os.system('rm *chunk*')

    # PARSER
    command = 'orthomclBlastParser blast.out adjusted_fasta >> similarSequences.txt'
    run_command(command)
    command = 'orthomclLoadBlast orthomcl.config similarSequences.txt'
    run_command(command)
    if os.path.isdir('pairs'):
        os.system('rm -r pairs/')
    command = 'orthomclPairs orthomcl.config sql.log cleanup=no'
    run_command(command)
    command = 'orthomclDumpPairsFiles orthomcl.config'
    run_command(command)

    # MCL
    command = 'mcl mclInput --abc -I 1.5 -o mclOutput > /dev/null'
    print command
    # by Sacha added bin_dir
    os.system(command)

    # GROUPS
    command = 'orthomclMclToGroups OG_ 1 < mclOutput > groups.txt'
    run_command(command)

    # SQL CLEANUP
    cursor = openDB(sql_user, sql_passwd)
    command = 'DROP SCHEMA %s;' % name
    print command
    cursor.execute(command)
    command = 'DROP USER %s;' % name
    print command
    cursor.execute(command)

    print 'Removing lockfile'
    os.system('rm %s' % lockname)
    cursor.close()

    print 'All set. Everything looks hunky-dory :) and the output in the groups.txt file in your working directory'
