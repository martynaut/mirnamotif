import click
from motiffind import motiffind, process_input


@click.command()
@click.argument('inputfile')
@click.option('--direction', '-d')
@click.option('--database', '-db')
@click.option('--results_folder', '-f')
@click.option('--negative_file', '-n')
def main(inputfile, direction, database, results_folder, negative_file):
    if not direction:
        direction = 'f'
    if not results_folder:
        results_folder = 'motif_search_results'
    if database:
        negative_file = './mirnamotif/neg_seq/miR_'+database+'.fa'
    if not database:
        database = 'loops'
    if not negative_file:
        negative_file = './mirnamotif/hairpin.fa'
    click.echo("Filename with sequences:, {}".format(inputfile))
    click.echo("Direction:, {}".format(direction))
    click.echo("Results folder: {}".format(results_folder))
    file = open(inputfile)
    data = file.read()
    list_of_seq = process_input(data, database)
    results_file, list_of_motifs, status = motiffind(list_of_seq=list_of_seq, direction=direction, folder=results_folder,
                                                     negative_file=negative_file)
    if status:
        click.echo("Results saved in {}".format(results_file))
        click.echo("Motifs found: {}".format(list_of_motifs))
    else:
        click.echo("No motifs were found")


if __name__ == "__main__":
    main()
