"""localizing loop of pre-miRNA from miRBase."""


import pandas as pd
import RNA


def find_arm(row):
    if row['-/+'] == '+':
        if row['start'] - row['start_pre'] < row['stop_pre'] - row['stop']:
            return '5p'
        else:
            return '3p'
    if row['-/+'] == '-':
        if row['start'] - row['start_pre'] < row['stop_pre'] - row['stop']:
            return '3p'
        else:
            return '5p'


def create_loc():
    coordinates_mmu = pd.read_csv('mmu.gff3', skiprows=13, sep='\t',
                                  names=['chr', 'col1', 'type', 'start', 'stop',
                                         'col2', '-/+', 'col3', 'ids'])
    coordinates_hsa = pd.read_csv('hsa.gff3', skiprows=13, sep='\t',
                                  names=['chr', 'col1', 'type', 'start', 'stop',
                                         'col2', '-/+', 'col3', 'ids'])
    coordinates_ath = pd.read_csv('ath.gff3', skiprows=13, sep='\t',
                                  names=['chr', 'col1', 'type', 'start', 'stop',
                                         'col2', '-/+', 'col3', 'ids'])
    coordinates = pd.concat([coordinates_mmu, coordinates_ath, coordinates_hsa])
    coordinates['ID'] = coordinates['ids'].str.extract('ID=(.*?);', expand=False)
    coordinates['Alias'] = coordinates['ids'].str.extract('Alias=(.*?);', expand=False)
    coordinates['Name'] = coordinates['ids'].str.extract('Name=(.*?)(?:$|;)', expand=False)
    coordinates['From'] = coordinates['ids'].str.extract('Derives_from=(.*)', expand=False)

    loc_info = coordinates
    sequences_ath = pd.read_csv('./miR_sh_ath', sep=' ', header=None, names=['Name', 'seq'])
    sequences_mmu = pd.read_csv('./miR_sh_mmu', sep=' ', header=None, names=['Name', 'seq'])
    sequences_hsu = pd.read_csv('./miR_sh_hsa', sep=' ', header=None, names=['Name', 'seq'])
    sequences = pd.concat([sequences_ath, sequences_hsu, sequences_mmu])
    pre_mirbase = loc_info[loc_info['type'] == 'miRNA_primary_transcript'].copy()

    matures = loc_info[loc_info['type'] == 'miRNA'].copy()
    pre_mirbase.drop(['col1', 'col2', 'col3', 'Alias', 'type', 'From', 'chr', 'ids'], inplace=True, axis=1)
    matures.drop(['col1', 'col2', 'col3', 'Alias', 'type', 'ids'], inplace=True, axis=1)
    matures = matures.join(pre_mirbase.set_index('ID'), on='From', how='left', rsuffix='_pre')
    matures['arm'] = matures['Name'].str[-2:]
    matures.loc[~matures['arm'].isin(['3p',
                                      '5p']), 'arm'] = \
        matures[~matures['arm'].isin(['3p',
                                      '5p'])].apply(lambda x: find_arm(x), axis=1)

    matures['Name'] = matures['Name_pre'].map(str) + matures['Name'].str.extract(
        '(\-[35]{1}p)', expand=False).fillna('')
    joined = pre_mirbase.join(matures.set_index('From'), on='ID', how='left', rsuffix='_mirna')
    localizations = pd.DataFrame(columns=['Name', 'seq'])

    joined = joined[~joined['arm'].isna()]
    joined = joined[joined['start'] <= joined['start_mirna']]
    joined = joined[joined['stop'] >= joined['stop_mirna']]
    temp = pd.DataFrame(joined.groupby('ID')['Name', 'stop_mirna', 'arm'].agg({'Name': 'count',
                                                                               'stop_mirna': 'count',
                                                                               'arm': 'sum'}))
    temp.columns = ['how many mature mirnas in pre', 'count', 'concat']
    both_mirnas = temp[(temp['concat'].str.contains('3p')) &
                       (temp['concat'].str.contains('5p'))]
    single_mirna = temp[~((temp['concat'].str.contains('3p')) &
                        (temp['concat'].str.contains('5p')))]
    joined_both = joined[joined['ID'].isin(both_mirnas.index.values)]
    joined_single = joined[joined['ID'].isin(single_mirna.index.values)]

    for premirna in joined_both['ID'].unique():

        data = joined_both[joined_both['ID'] == premirna].reset_index()
        sequence = sequences[sequences['Name'] == data['Name'][0]]['seq'].values[0]
        if data.shape[0] > 2:
            data_5_plus = data[(data['arm'] == '5p') & (data['-/+'] == '+')]
            data_5_minus = data[(data['arm'] == '5p') & (data['-/+'] == '-')]
            data_3_plus = data[(data['arm'] == '3p') & (data['-/+'] == '+')]
            data_3_minus = data[(data['arm'] == '3p') & (data['-/+'] == '-')]
            data_5_plus = data_5_plus.groupby('arm').agg({'start': 'first',
                                                          'stop': 'first',
                                                          '-/+': 'first',
                                                          'ID': 'first',
                                                          'Name': 'first',
                                                          'chr': 'first',
                                                          'start_mirna': 'min',
                                                          'stop_mirna': 'max',
                                                          })

            data_5_minus = data_5_minus.groupby('arm').agg({'start': 'first',
                                                            'stop': 'first',
                                                            '-/+': 'first',
                                                            'ID': 'first',
                                                            'Name': 'first',
                                                            'chr': 'first',
                                                            'start_mirna': 'min',
                                                            'stop_mirna': 'max',
                                                            })

            data_3_plus = data_3_plus.groupby('arm').agg({'start': 'first',
                                                          'stop': 'first',
                                                          '-/+': 'first',
                                                          'ID': 'first',
                                                          'Name': 'first',
                                                          'chr': 'first',
                                                          'start_mirna': 'min',
                                                          'stop_mirna': 'max',
                                                          })

            data_3_minus = data_3_minus.groupby('arm').agg({'start': 'first',
                                                            'stop': 'first',
                                                            '-/+': 'first',
                                                            'ID': 'first',
                                                            'Name': 'first',
                                                            'chr': 'first',
                                                            'start_mirna': 'min',
                                                            'stop_mirna': 'max',
                                                            })
            data = pd.concat([data_3_minus, data_3_plus, data_5_minus, data_5_plus]).reset_index()

        coordinates = []
        coordinates += list(data['start'].unique())
        coordinates += list(data['stop'].unique())
        coordinates += list(data['start_mirna'].unique())
        coordinates += list(data['stop_mirna'].unique())
        coordinates = list(map(int, coordinates))
        coordinates = sorted(coordinates)

        if data.loc[0, '-/+'] == '+':
            localizations = localizations.append(
                {'Name': data.loc[0, 'Name'],
                 'seq': sequence[coordinates[2] - coordinates[0] + 1: coordinates[3] - coordinates[0]]
                 }, ignore_index=True
            )

        elif data.loc[0, '-/+'] == '-':

            localizations = localizations.append(
                {'Name': data.loc[0, 'Name'],
                 'seq': sequence[coordinates[-1] - coordinates[-3] + 1: coordinates[-1] - coordinates[-4]]
                 }, ignore_index=True
            )

        else:
            print('not all records has the same orientation')
            return 1

    for premirna in joined_single['ID'].unique():
        temp_coordinates = []
        data = joined_single[joined_single['ID'] == premirna].reset_index()
        sequence = sequences[sequences['Name'] == data['Name'][0]]['seq'].values[0]
        data = data.groupby('arm').agg({'start': 'first',
                                        'stop': 'first',
                                        '-/+': 'first',
                                        'ID': 'first',
                                        'Name': 'first',
                                        'chr': 'first',
                                        'start_mirna': 'min',
                                        'stop_mirna': 'max',
                                        }).reset_index()

        coordinates = []
        coordinates += list(data['start'].unique())
        coordinates += list(data['stop'].unique())
        coordinates += list(data['start_mirna'].unique())
        coordinates += list(data['stop_mirna'].unique())
        coordinates = list(map(int, coordinates))
        coordinates = sorted(coordinates)

        if (abs(coordinates[0] - coordinates[1]) < abs(coordinates[2] - coordinates[3])
                and data.loc[0, '-/+'] == '+'):
            known = 5
        elif (abs(coordinates[0] - coordinates[1]) < abs(coordinates[2] - coordinates[3])
              and data.loc[0, '-/+'] == '-'):
            known = 3
        elif data.loc[0, '-/+'] == '+':
            known = 3
        else:
            known = 5

        structure = RNA.fold(sequences[sequences['Name'] == data['Name'][0]]['seq'].values[0])[0]

        if data.loc[0, '-/+'] == '+':
            if known == 5:
                predicted_start = coordinates[1] - coordinates[0] + 1
                pre_mirna = str(structure[:predicted_start]).count('(')
                mirna_len = coordinates[2] - coordinates[1] + 1
                mirna_len_com = str(
                    structure[predicted_start:predicted_start+mirna_len+1]).count('(')
                closing = 0
                start_second_strand = coordinates[-1] - coordinates[0] + 2
                while closing < pre_mirna:
                    start_second_strand = start_second_strand - 1
                    closing = str(structure[start_second_strand:]).count(')')
                closing2 = 0
                start_second_strand2 = coordinates[-1] - coordinates[0] + 2
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 - 1
                    closing2 = str(structure[start_second_strand2:]).count(')')

                stop_second_strand = start_second_strand2
                start_second_strand = start_second_strand + 1 + coordinates[0]
                stop_second_strand = stop_second_strand + 2 + coordinates[0]
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()

            if known == 3:
                predicted_start = coordinates[2] - coordinates[0] + 1
                pre_mirna = str(structure[predicted_start:]).count(')')

                mirna_len = coordinates[2] - coordinates[1] + 1
                mirna_len_com = str(
                    structure[predicted_start-mirna_len:predicted_start]).count(')')

                closing = 0
                start_second_strand = 0
                while closing < pre_mirna:
                    start_second_strand = start_second_strand + 1
                    closing = str(structure[:start_second_strand]).count('(')
                closing2 = 0
                start_second_strand2 = 0
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 + 1
                    closing2 = str(structure[:start_second_strand2]).count('(')
                stop_second_strand = start_second_strand2
                start_second_strand = start_second_strand + 1 + coordinates[0]
                stop_second_strand = stop_second_strand + 1 + coordinates[0]
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()

        if data.loc[0, '-/+'] == '-':
            if known == 5:
                predicted_start = coordinates[3] - coordinates[2]
                pre_mirna = str(structure[:predicted_start]).count('(')
                mirna_len = coordinates[2] - coordinates[1] - 1
                mirna_len_com = str(
                    structure[predicted_start:predicted_start+mirna_len+1]).count('(')

                closing = 0
                start_second_strand = coordinates[3] - coordinates[0] + 1
                while closing < pre_mirna:
                    start_second_strand = start_second_strand - 1
                    closing = str(structure[start_second_strand:]).count(')')
                closing2 = 0
                start_second_strand2 = coordinates[3] - coordinates[0] + 1
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 - 1
                    closing2 = str(structure[start_second_strand2:]).count(')')

                stop_second_strand = start_second_strand2
                start_second_strand = coordinates[3] - start_second_strand - 1
                stop_second_strand = coordinates[3] - stop_second_strand - 1
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()

            if known == 3:
                predicted_start = coordinates[3] - coordinates[1]
                pre_mirna = str(structure[predicted_start:]).count(')')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure[predicted_start-mirna_len:predicted_start]).count(')')

                closing = 0
                start_second_strand = 0
                while closing < pre_mirna+1:
                    start_second_strand = start_second_strand + 1
                    closing = str(structure[:start_second_strand]).count('(')
                if pre_mirna == 0:
                    start_second_strand = 0
                closing2 = 0
                start_second_strand2 = 0
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 + 1
                    closing2 = str(structure[:start_second_strand2]).count('(')

                stop_second_strand = start_second_strand2
                start_second_strand = coordinates[3] - start_second_strand - 1
                stop_second_strand = coordinates[3] - stop_second_strand - 1
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()

        coordinates = []
        coordinates += list(data['start'].unique())
        coordinates += list(data['stop'].unique())
        coordinates += list(data['start_mirna'].unique())
        coordinates += list(data['stop_mirna'].unique())
        coordinates += temp_coordinates
        coordinates = list(map(int, coordinates))
        coordinates = sorted(coordinates)

        if data.loc[0, '-/+'] == '+':

            localizations = localizations.append(
                {'Name': data.loc[0, 'Name'],
                 'seq': sequence[coordinates[2] - coordinates[0] + 1: coordinates[3] - coordinates[0]]
                 }, ignore_index=True
            )

        elif data.loc[0, '-/+'] == '-':
            localizations = localizations.append(
                {'Name': data.loc[0, 'Name'],
                 'seq': sequence[coordinates[-1] - coordinates[-3] + 1: coordinates[-1] - coordinates[-4]]
                 }, ignore_index=True
            )
        else:
            print('not all records has the same orientation')
            return 1
    localizations = localizations.sort_values('Name')
    localizations[localizations['Name'].str.startswith('hsa')].to_csv('./miR_linking_hsa', sep=' ',
                                                                      index=False, header=False)
    localizations[localizations['Name'].str.startswith('mmu')].to_csv('./miR_linking_mmu', sep=' ',
                                                                      index=False, header=False)
    localizations[localizations['Name'].str.startswith('ath')].to_csv('./miR_linking_ath', sep=' ',
                                                                      index=False, header=False)
