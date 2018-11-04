import click
from motifsearch import search


@click.command()
@click.argument('seq')
@click.option('--seq2', '-s')
@click.option('--seq_type', '-t')
@click.option('--direction', '-d')
@click.option('--overlap', '-o')
@click.option('--database', '-db')
@click.option('--file_prefix', '-fp')
def main(seq, seq2='', seq_type='loops', direction='f', overlap=True, database='hsa', file_prefix='local'):
    if not seq2:
        seq2 = ''
    if not seq_type:
        seq_type = 'loops'
    if not direction:
        direction = 'f'
    if not overlap or overlap == 'True':
        overlap = True
    else:
        overlap = False
    if not database:
        database = 'hsa'
    if not file_prefix:
        file_prefix = 'local'

    if seq2:
        click.echo("Motifs searched: {}, {}".format(seq, seq2))
    else:
        click.echo("Motif searched: {}".format(seq))
    click.echo("Type of sequences searched: {}".format(seq_type))
    click.echo("Direction:, {}".format(direction))
    click.echo("Overlap: {}".format(overlap))
    click.echo("Database: {}".format(database))
    click.echo("sequences searched: {}, {}".format(seq, seq2))
    filenames, status = search(seq1=seq, seq2=seq2, seq_type=seq_type, direction=direction,
       overlap=overlap, database=database, email=file_prefix)
    if status:
        click.echo("Results saved in {}\nand {}".format(filenames[0], filenames[0]))
    else:
        click.echo("No matching sequences found")


if __name__ == "__main__":
    main()
