# Standard library imports
import io
from statistics import stdev

# Third-party imports
from flask import Flask, request, send_file, render_template as rt, Response
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import psycopg2

# Local imports
from cog import Cog
from protein import Protein

# Global variables
app = Flask(__name__)
conn = None


def get_db() -> None:
    """
    The function checks if there is a connection saved in the global variable
    'conn'. If not, it will create a connection and save it in the global
    variable 'conn'.
    """
    global conn
    if conn is None:
        conn = psycopg2.connect(
            host='145.97.18.224',
            user='bpcoogsa1',
            password='bpcoogsa1',
            database='bpcoogsa1_db')


def get_organisms() -> dict:
    """
    The function gets a database connection and queries the database for
    the organism_id and latin_name from the organisms table.

    :return: dictionary with organism_id as key and latin_name as value.
    """
    get_db()
    cur = conn.cursor()
    cur.execute(f'SELECT organisme_id, latijnse_naam FROM organismen')
    organisms = cur.fetchall()
    return {org_id: latin_name for org_id, latin_name in organisms}


def fetch_cogs(cog_filter) -> list[Cog]:
    """
    The function gets a database connection and queries the database for
    'cog_id', 'go_annotation' and 'function' from the table 'cogs' and then for
    'protein_id', 'organism_id', 'sequence' and 'cog_id' from the joined tables
    of 'proteins' and 'proteins_cogs'.

    A dictionary of the latin names of all organisms in the database is fetched
    and used to create a Protein object for every protein in the database.
    Subsequently, the Protein objects are used to create a Cog object for every
    COG in the database. The apply_filters() function is called if the
    cog_filter argument has been filled.

    :param cog_filter: tuple containing the three ranges for the COG filter.
    :return: list of Cog objects.
    """
    # Get database connection and cursor.
    get_db()
    cur = conn.cursor()

    # Query database for COG and protein data.
    cur.execute(f'SELECT cog_id, go_annotation, functie FROM cogs;')
    cog_data = cur.fetchall()
    cur.execute(f'SELECT eiwitten.eiwit_id, organisme_id, sequentie, cog_id '
                f'FROM eiwitten JOIN eiwitten_cogs '
                f'ON eiwitten.eiwit_id = eiwitten_cogs.eiwit_id;')
    protein_data = cur.fetchall()

    # Get the organism dictionary and create Protein and Cog objects.
    org_dict = get_organisms()
    proteins = [Protein(p[0], org_dict.get(p[1]), p[3], p[2])
                for p in protein_data]
    cogs = [Cog(cd[0], cd[1], cd[2],
                [p for p in proteins if p.get_cog_id() == cd[0]])
            for cd in cog_data]

    # Return filtered or unfiltered COGs depending on the cog_filter argument.
    if cog_filter:
        return apply_filters(cogs, cog_filter[0], cog_filter[1], cog_filter[2])
    else:
        return cogs


def get_single_cog(cog_id) -> list[tuple]:
    """
    The function gets a database connection and queries the database for
    'cog_id', 'go_annotation' and 'function' from the table 'cogs' and then for
    'protein_id', 'organism_id', 'sequence' and 'cog_id' from the joined tables
    of 'proteins' and 'proteins_cogs' and fetches the data where the 'cog_id'
    corresponds to the argument given to the function.

    :param cog_id: string identifier for the COG to fetch.
    :return: list of tuples with the data for each protein in the COG.
    """
    get_db()
    cur = conn.cursor()
    cur.execute(f"SELECT eiwitten.eiwit_id, organisme_id, sequentie, cog_id "
                f"FROM eiwitten JOIN eiwitten_cogs "
                f"ON eiwitten.eiwit_id = eiwitten_cogs.eiwit_id "
                f"WHERE cog_id LIKE '{cog_id}'")
    protein_data = cur.fetchall()
    cur.execute(f"SELECT go_annotation FROM cogs "
                f"WHERE cog_id LIKE '{cog_id}'")
    go_annotation = cur.fetchall()[0][0]
    return protein_data, go_annotation


def plot_cogs(cogs) -> str:
    """
    The function creates a matplotlib figure containing a scatterplot of the
    Cog objects in the cogs list. The figure is then converted for display in
    the web application.

    :param cogs: list of Cog objects to display in the plot.
    :return: string with the HTML code for the image.
    """
    # Call zip function to prepare data for plotting.
    x_vals, y_vals, sizes = zip(*[(float(cog.get_avg_seq_len()),
                                   float(cog.get_rsd()),
                                   int(cog.get_methionines().split('/')[0][1]))
                                  for cog in cogs])

    # Create a Figure and FigureCanvas and configure the plot.
    fig = Figure()
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.scatter(x_vals, y_vals, sizes=sizes, alpha=0.5, c="#5271ff")
    ax.set_ylabel('RSD of average sequence length')
    ax.set_xlabel('Average sequence length')
    ax.set_title(f'{len(cogs)} COGs in database')
    ax.grid(True)

    # Convert the plot for display in web application.
    img = io.StringIO()
    fig.savefig(img, format='svg')
    return '<svg' + img.getvalue().split('<svg')[1]


def apply_filters(cogs, rsd, seq_len, met) -> list[Cog]:
    """
    The function recreates the list of Cog objects using a list comprehension
    for each filter argument that is filled. The resulting list of Cog objects
    is returned.

    :param cogs: list of Cog objects.
    :param rsd: tuple containing the upper and lower bounds for the relative
    standard deviation.
    :param seq_len: tuple containing the upper and lower bounds for the average
    sequence length.
    :param met: tuple containing the upper and lower bounds for the number of
    methionines as first amino acid.
    :return: list of Cog objects.
    """
    if rsd:
        cogs = [cog for cog in cogs if float(rsd[0])
                <= float(cog.get_rsd()) <= float(rsd[1])]
    if seq_len:
        cogs = [cog for cog in cogs if float(seq_len[0])
                <= float(cog.get_avg_seq_len()) <= float(seq_len[1])]
    if met:
        cogs = [cog for cog in cogs if int(met[0])
                <= int(cog.get_methionines().split('/')[0][1]) <= int(met[1])]
    return cogs


def barchart_proteins(proteins) -> str:
    """
    The function creates a matplotlib figure containing a barchart of the
    Protein objects in the proteins list. The figure is then converted for
    display in the web application.

    :param proteins: list of Protein objects to display in the plot.
    :return: string with the HTML code for the image.
    """
    # Parse proteins for data to display in the barchart.
    protein_ids = [p.get_protein_id() for p in proteins]
    sequence_lengths = [len(p.get_sequence()) for p in proteins]
    colour_dict = {0: 'tab:orange', 1: '#5271FF'}
    bar_colors = [colour_dict.get(int(p.get_sequence()[0] == 'M'))
                  for p in proteins]

    # Create a Figure and FigureCanvas and set up the data and barchart.
    fig, ax = plt.subplots()
    FigureCanvas(fig)
    ax.bar(protein_ids, sequence_lengths, color=bar_colors)
    ax.set_ylabel('Sequence length')
    ax.set_title('Proteins in COG')
    ax.legend(title='Protein id')
    plt.xticks(rotation=90)
    plt.gcf().subplots_adjust(bottom=0.35)

    # Convert the plot for display in web application.
    img = io.StringIO()
    fig.savefig(img, format='svg')
    return '<svg' + img.getvalue().split('<svg')[1]


def write_fasta(cog_number, proteins) -> None:
    """
    The function opens a .fasta file using cog_number as the name to write the
    FASTA sequence for each protein in the proteins list to it. This is done
    using the __repr__() function that is overwritten in the Protein class.

    :param cog_number: string identifier for the COG.
    :param proteins: list of Protein objects.
    """
    with open(f'{cog_number}.fasta', 'w') as output:
        for p in proteins:
            output.write(p.__repr__())


@app.route('/', methods=['GET'])
def login() -> str:
    """
    The function renders the login page where users can connect to a database.
    """
    return rt('login.html')


@app.route('/help', methods=['GET'])
def help_page() -> str:
    """
    The function renders the help page where users can learn about COGselect.
    """
    return rt('help.html')


@app.route('/db_overview', methods=['POST'])
def display_data() -> str:
    """
    The function fetches Cog objects and finds the ranges for their relative
    standard deviation, methionines as first amino acid and average sequence
    length. The Cog objects are displayed in a scrollbox, plotted and their
    values are used to fill the filter form.
    """
    # Fetch Cog objects and store the data used to filter them.
    cogs = fetch_cogs(cog_filter=False)
    methionines = [cog.get_methionines().split('/')[0][1] for cog in cogs]
    rsd = [cog.get_rsd() for cog in cogs]
    avg_len = [cog.get_avg_seq_len() for cog in cogs]
    return rt('db_overview.html',
              cogs=cogs,
              plot=plot_cogs(cogs),
              rsd_range=f'{min(rsd)} - {max(rsd)}',
              meth_range=f'{min(methionines)} - {max(methionines)}',
              avg_len_range=f'{min(avg_len)} - {max(avg_len)}')


@app.route('/filter_cogs', methods=['POST'])
def display_filtered() -> str:
    """
    The function renders the database overview page with a filtered list of Cog
    objects. The filter values used are loaded back up to the filter form.

    :return:
    """
    # Get the values used as filters and format them into a tuple of the upper
    # and lower bounds for each range.
    rsd = request.form['rsd'].strip(' ').split('-')
    avg_length = request.form['avg_len'].strip(' ').split('-')
    methionines = request.form['meth_starts'].strip(' ').split('-')
    cog_filter = (tuple(boundary for boundary in rsd),
                  tuple(boundary for boundary in avg_length),
                  tuple(boundary for boundary in methionines))

    # Cog objects are fetched using the cog_filter and the database overview is
    # rendered.
    cogs = fetch_cogs(cog_filter)
    return rt('db_overview.html',
              title='COGs in database',
              cogs=cogs,
              plot=plot_cogs(cogs),
              rsd_range=f'{cog_filter[0][0]} - {cog_filter[0][1]}',
              avg_len_range=f'{cog_filter[1][0]} - {cog_filter[1][1]}',
              meth_range=f'{cog_filter[2][0]} - {cog_filter[2][1]}')


@app.route('/view_cog', methods=['GET', 'POST'])
def view_cog() -> str:
    """
    The function gets the cog_id from the request form and fetches data using
    the cog_id. Protein objects are created and analyzed for display on the
    webpage.

    :return: string that renders the webpage.
    """
    cog_id = request.form['cog_id']
    protein_data, go_annotation = get_single_cog(cog_id)
    go_link = f'http://amigo.geneontology.org/amigo/term/GO:{go_annotation}'
    org_dict = get_organisms()
    proteins = [Protein(p[0], org_dict.get(p[1]), p[3], p[2])
                for p in protein_data]
    protein_lengths = [p.get_length() for p in proteins]
    write_fasta(cog_id, proteins)
    return rt('cog.html',
              cog_id=cog_id,
              go_link=go_link,
              proteins=proteins,
              plot=barchart_proteins(proteins),
              avg_seq_length=round(sum(protein_lengths) / len(proteins), 2),
              std_avg_seq_length=round(float(stdev(protein_lengths)), 2))


@app.route('/download', methods=['POST'])
def download_fasta() -> Response:
    """
    The function gets the cog_id from the request form and uses it to generate
    a file path to return using Flasks send_file. The file is then downloaded.

    :return: Flask Response to download file.
    """
    cog_id = request.form['cog_id']
    path = f'{cog_id}.fasta'
    return send_file(path, as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
