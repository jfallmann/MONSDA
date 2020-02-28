#!/usr/bin/python3

### execute the script applying: python tRNAdb.py > ../tRNAdb/tRNAdb.fasta
### Instruction were found under: http://zetcode.com/db/mysqlpython/

### DESCRIPTION: Retrieving data from the database

#import MySQLdb as mdb
import pymysql
#import PyMySQL
#print("Here")

#import sys

# connecting to the database with the parameters: server, user, password, db
#con = PyMySQL.connect("maitai.bioinf.uni-leipzig.de", "anneh", "HyhAttDB11","trna")
db = pymysql.connect("maitai.bioinf.uni-leipzig.de", "anneh", "HyhAttDB11","trna")

cursor = db.cursor()
cursor.execute("SELECT * FROM rna")

rows = cursor.fetchall()
for row in rows:
    print (">","t_",row[0],row[1],row[2],row[3],row[4],"\n",row[7])

cursor.execute("SELECT * FROM gdb")

rows = cursor.fetchall()
for row in rows:
    print (">","g_",row[0],row[1],row[2],row[3],row[4],"\n",row[7])

db.close()
