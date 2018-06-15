# Version:  1.0
# Date:     03-06-2018
# Author:   Project group 10
# Function: This module sets up a connection to the mysql database with the connect function.
#           Furthermore, it accepts BLAST results in the form of read and hit data in the insert_read_and_data function,
#           which is then passed along to the insert_read and insert_hits functions respectively to insert these
#           results into the database.

import mysql.connector as mscon


def connect(username, hostname, databasename, password):
    """This function makes a connection with the mysql database.

    This function accepts a username, hostname, databasename and password and tries to connect with
    the mysql database. If successfull, it returns the connection object and a cursor object.
    If an error occurs, it prints what error occurred and returns False.
    """
    try:
        con = mscon.connect(user=username, host=hostname, database=databasename, password=password)
    except mscon.ProgrammingError as e:
        print("Could not connect to the database, the following error occurred:\n{}".format(e))
        return False
    else:
        cursor = con.cursor()
    return con, cursor


def insert_read_and_data(cursor, read_list, hit_list):
    """Accepts a cursor, read data and hit data and passes it to insert_read and insert_hits.

    This function accepts a cursor, a read list with header, sequence and quality score and a 2d list with
    data about the hits. It checks if the read is already in the database with the insert_read function, which inserts
    it if it isn't, and then inserts the corresponding hits with the insert_hits function if the read wasn't already
    in the database.
    """
    read_id = insert_read(cursor, read_list)
    if read_id:
        insert_hits(cursor, read_id, hit_list)
    else:
        print("Insertion failed.")


def insert_read(cursor, read_list):
    """Accepts a cursor and list with read info, inserts the read into the database if it doesnt exist and returns its
    DNA_READ_id. If it already existed, returns False

    This function accepts a cursor and a read list with header, sequence and quality score. If the header is already in
    the database, it returns false. If not, it inserts the read into the database. It then checks if there is a
    corresponding forward or reverse read already in the database. If there is, it updates the For_rev_ids of both reads
    so they reference each other.
    """
    new_header = read_list[0]
    id_query = """
      SELECT DNA_READ_id FROM DNA_READ
      WHERE Header = '{}'
    """
    cursor.execute(id_query.format(new_header))
    if not cursor.fetchall():
        read_query = """
          Insert INTO DNA_READ (Header, Sequentie, Quality_score)
          Values('{}',
                 '{}',
                 {}
        );
        """.format(new_header, read_list[1], read_list[2])

        cursor.execute(read_query)
        cursor.execute(id_query.format(new_header))
        new_entry_id = cursor.fetchall()[0][0]

        if new_header[-1] == "1":
            inverse_header = new_header[0:-1]+"2"
        elif new_header[-1] == "2":
            inverse_header = new_header[0:-1] + "1"
        else:
            inverse_header = False

        if inverse_header:
            cursor.execute(id_query.format(inverse_header))
            result_list = cursor.fetchall()
            if result_list:
                old_entry_id = result_list[0][0]
                read_update_query = """
                  UPDATE DNA_READ
                  SET For_rev_id = {}
                  WHERE Header = '{}'"""
                cursor.execute(read_update_query.format(old_entry_id, new_header))
                cursor.execute(read_update_query.format(new_entry_id, inverse_header))
        return new_entry_id
    else:
        print("That header already exists")
        return False


def insert_hits(cursor, read_id, hit_list):
    """Accepts a cursor, a read_id and 2d list with hit data and inserts the hit data into the database.

    This function accepts a cursor, a read_id and 2d list with hit data, each sublist containing a BLAST score,
    percentage query cover, percentage identity, percentage positives, E value, organism name, protein name,
    protein comment and protein access code. For each hit, it checks if the organism and protein are already in
    the database and inserts them if they're not, after which their ID's are retrieved. The remaining hit data is then
    inserted into the hit table.
    """
    exists_query = """
      SELECT {} FROM {}
      WHERE {} = '{}'"""

    for hit in hit_list:
        organisme = hit[5].replace("\'", "\'\'")
        cursor.execute(exists_query.format('ORGANISME_id', 'ORGANISME', 'Organisme_naam', organisme))
        results = cursor.fetchall()
        if not results:
            cursor.execute("""
                             INSERT INTO ORGANISME (Organisme_naam)
                             VALUES ('{}')""".format(organisme))
            cursor.execute(exists_query.format('ORGANISME_id', 'ORGANISME', 'Organisme_naam', organisme))
            results = cursor.fetchall()
        organisme_id = results[0][0]

        access_code = hit[8]
        cursor.execute(exists_query.format('EIWIT_id', 'EIWIT', 'Accessiecode', access_code))
        results = cursor.fetchall()
        if not results:
            if hit[7] != "":
                cursor.execute("""
                  INSERT INTO EIWIT (Naam, Eiwit_comment, Accessiecode)
                  VALUES ('{}', '{}', '{}')""".format(hit[6].replace("\'", "\'\'"), hit[7].replace("\'", "\'\'"), access_code))
            else:
                cursor.execute("""
                  INSERT INTO EIWIT (Naam, Eiwit_comment, Accessiecode)
                  VALUES ('{}', Null, '{}')""".format(hit[6].replace("\'", "\'\'"), access_code))
            cursor.execute(exists_query.format('EIWIT_id', 'EIWIT', 'Accessiecode', access_code))
            results = cursor.fetchall()
        eiwit_id = results[0][0]

        cursor.execute("""
          INSERT INTO HIT (DNA_READ_id, ORGANISME_id, EIWIT_id, Score, Query_cover, Identity, Positives, E_value)
          VALUES ({}, {}, {}, {}, {}, {}, {}, {})""".format(read_id,
                                                            organisme_id,
                                                            eiwit_id,
                                                            hit[0],
                                                            hit[1],
                                                            hit[2],
                                                            hit[3],
                                                            hit[4]))
