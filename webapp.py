# Version:  1.0
# Date:     14-06-2018
# Author:   Project group 10
# Function: This module connects with the mysql database in which the BLAST data is stored
#           and visualises it in a web application. It has a homepage in which a graph of the most
#           found organisms is shown and has the options to display a graph with the most found proteins
#           and to search the database and display the results in an organised manner.

from flask import Flask, render_template, request, redirect, url_for
import mysql.connector as mscon
import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from io import BytesIO
import base64

app = Flask(__name__)

# Dictionary with info about the graph
GRAPHS = { 'soorten':
{       'html': '',
        'title': 'Soorten',
        'description': 'Meest voorkomende soorten'},
          'eiwitten':
 {      'html': '',
        'title': 'Soorten',
        'description': 'Meest voorkomende soorten'}
          }

def makecon():
    # makes connection with mySQL database
    #returns cursor
    con = mscon.connect(user="sql7241825", host="sql7.freemysqlhosting.net", database="sql7241825", password="7cCBxT27sc")
    cursor = con.cursor()
    return cursor, con

def getdata(cursor):
    #input: cursor
    #gets organism and protein data from database
    #returns list with data

    table = 'ORGANISME'
    data = {}

    query = "select * from {}".format(table)
    cursor.execute(query)
    sub = []
    for value in cursor:
        sub.append(value[-1])
    data[table] = sub

    return data

def makegraph(cursor):
    #input: cursor
    #creates barcharts of most common organisms
    #returns list of most common organisms
    querys = [
        "SELECT count(*), Organisme_naam FROM HIT natural join ORGANISME group by ORGANISME_id order by count(*) desc",
        "SELECT count(*), naam FROM HIT natural join EIWIT group by Naam order by count(*) desc"]
    mostcommon = {
        'soorten': [],
        'eiwitten': []
    }
    for i in range(2):
        cursor.execute(querys[i])
        x = []
        y = []
        for value in cursor:
            x.append(value[1])
            y.append(value[0])
        if i == 0:
            id = 'soorten'
            for index in range(10):
                mostcommon[id].append(x[index])
        if i == 1:
            id = 'eiwitten'
            for index in range(10):
                mostcommon[id].append(x[index])

        plt.bar(range(10), y[:10], color='#5f021f', align='center')
        plt.xticks(range(10), range(1, 11))
        plt.title("Meest voorkomende " + id)
        plt.xlabel("")
        plt.ylabel("Frequentie")
        plt.tight_layout()
        graphfile = BytesIO()
        plt.savefig(graphfile, format='png')
        graphfile.seek(0)
        figdata_png = graphfile.getvalue()
        figdata_png = base64.b64encode(figdata_png)
        figdata_png = figdata_png.decode('utf-8')
        GRAPHS[id]['html'] = figdata_png
        plt.close()
    return mostcommon


def result_retriever(organism, protein, protein_comment, read_quality):
    """This function accepts an organism name, protein name, protein comment and read quality. It then queries
    the mysql database to find all results in which the organism name, protein name and protein comment are present
    and where the minimal read quality was the read quality given to the function. It then returns a list of tuples
    containing the data of each result.
    """
    cursor, con = makecon()

    if read_quality == "":
        read_quality = 0
    else:
        try:
            read_quality = int(read_quality)
        except ValueError:
            read_quality = 0

    organism = '%'+organism+'%'
    protein = '%'+protein+'%'
    protein_comment ='%'+protein_comment+'%'


    query = """
      SELECT h.HIT_id, o.Organisme_naam, e.Naam, h.Score, h.Query_cover, h.Identity, h.E_value FROM HIT h
      NATURAL JOIN ORGANISME o
      NATURAL JOIN EIWIT e
      NATURAL JOIN DNA_READ r
      WHERE o.Organisme_naam LIKE %s
      AND e.Naam LIKE %s
      AND e.Eiwit_comment LIKE %s
      AND r.Quality_score > %s
      ORDER BY h.Score desc, h.Query_cover desc, h.Identity desc"""
    cursor.execute(query, (organism, protein, protein_comment, read_quality))
    result_list = []

    for element in cursor:
        try:
            protein_name = element[2].split('[')[0]
            result_list.append((element[0], element[1], protein_name, element[3], element[4], element[5], element[6]))
        except IndexError:
            result_list.append(element)
    cursor.close()
    con.close()
    return result_list


def get_hit_data(hit_id):
    """This function accepts a hit id and then queries the mysql database to retrieve all relevant information
    corresponding to that hit. It returns a 2d hit info list containing labels and values for all the data retrieved.
    It also returns a list of tuples containing the data of the read that the hit was found with, including the other
    hits that were found with that read. A similar list is returned for the reverse or forward read corresponding to
    that read.
    """
    cursor, con = makecon()
    query = """
      SELECT r.DNA_READ_id, r.For_rev_id, e.Naam, o.Organisme_naam, h.Score, h.Query_cover, h.Identity, h.Positives, h.E_value,
      e.Eiwit_comment, e.Accessiecode
      FROM DNA_READ r
      NATURAL JOIN HIT h
      NATURAL JOIN ORGANISME o
      NATURAL JOIN EIWIT e
      WHERE h.HIT_id = %s"""
    cursor.execute(query, (hit_id,))

    for element in cursor:
        result_tuple = element


    query = """
    SELECT r.Header, r.Quality_score, e.Naam, h.HIT_id, o.Organisme_naam FROM DNA_READ r
    NATURAL LEFT OUTER JOIN HIT h
    NATURAL LEFT OUTER JOIN EIWIT e
    NATURAL LEFT OUTER JOIN ORGANISME o
    WHERE r.DNA_READ_id = %s
    ORDER BY e.Naam, o.Organisme_naam"""
    cursor.execute(query, (result_tuple[0],))

    read_list1 = []
    for element in cursor:
        read_list1.append(element)

    read_list2 = []

    if result_tuple[1] is not None:

        cursor.execute(query, (result_tuple[1],))
        read_list2 = []
        for element in cursor:
            read_list2.append(element)

    label_list = ["Eiwit naam: ", "Organisme naam: ", "Bitscore: ", "Percentage query coverage: ",
                  "Percentage identity: ", "Percentage positives: ", "E value: ", "Eiwit comment: ", "Accessiecode: "]
    hit_info_list = []
    for value, label in zip(result_tuple[2:], label_list):
        hit_info_list.append([label, value])
    cursor.close()
    con.close()

    return hit_info_list, read_list1, read_list2


@app.route('/')
def gotohomepage():
    """This function redirects to the homepage."""
    return redirect(url_for('barchart',id='soorten'))

@app.route('/<id>')
def barchart(id):
    """This function renders the HTML template for the homepage of the application, showing a bar chart
    with the most found species or proteins, depending on the id that was given in the url."""
    cursor, con = makecon()
    graph = GRAPHS.get(id)
    mostcommon = makegraph(cursor)
    most10 = mostcommon.get(id)
    cursor.close()
    con.close()
    return render_template("graph.html", graph=graph, mostcommon=most10, id=id)


@app.route('/results', methods=['POST', 'GET'])
def results():
    """This function calls the result_retriever function to get results out of the mysql database
    based on the search information that was passed on from an HTML form. It then renders the HTML template
    that displays these results."""
    if request.method == 'POST':
        organism = request.form['organism']
        protein = request.form['protein']
        comment = request.form['comment']
        read_quality = request.form['read_quality']
        result_list = result_retriever(organism, protein, comment, read_quality)
        return render_template('results.html', result_list=result_list)
    else:
        return render_template('results.html', result_list=[[]])


@app.route('/hit/<hit_id>')
def hit(hit_id):
    """This function calls the get_hit_data function to retrieve all relevant hit information based on a hit id
    passed on in the URL. It the renders the HTML template that displays this hit information."""
    hit_info_list, read_list1, read_list2 = get_hit_data(hit_id)
    return render_template("hit_info.html", hit_info_list=hit_info_list, read_list1=read_list1, read_list2=read_list2)


if __name__ == '__main__':
    app.run()
