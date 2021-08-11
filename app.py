import base64
import os
from io import StringIO

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from CGRtools.exceptions import InvalidAromaticRing
from CGRtools.files import SDFRead, SMILESRead, SDFWrite
from stqdm import stqdm


def standardization_step(molecule, index, std_log, mistakes_log):
    molecule.clean2d()
    tmp_molecule = molecule.copy()
    standardized = False
    try:
        k = tmp_molecule.kekule()
        if tmp_molecule.check_valence():
            mistakes_log['indeces'].append(index)
            mistakes_log['mistakes'].append('valence_error')
            tmp_molecule.meta['standardized'] = standardized
            tmp_molecule.meta['mistake'] = True
            return molecule, std_log, mistakes_log
        else:
            s = tmp_molecule.standardize(fix_stereo=False)
            std_log['functional group'] += int(s)
            h = tmp_molecule.implicify_hydrogens(fix_stereo=False)
            std_log['implicify hydrogens'] += int(h)
            t = tmp_molecule.thiele()
            std_log['aromatization'] += int(t and not k)
            # c = tmp_molecule.standardize_charges(prepare_molecule=False)
            # std_log['charges'] += int(c)
            standardized = s or h or (t and not k)  # or c
    except InvalidAromaticRing:
        mistakes_log['indeces'].append(index)
        mistakes_log['mistakes'].append('invalid_aromatic_ring')
        tmp_molecule.meta['standardized'] = standardized
        tmp_molecule.meta['mistake'] = True
        return molecule, std_log, mistakes_log

    tmp_molecule.meta['standardized'] = standardized
    tmp_molecule.meta['mistake'] = False
    return tmp_molecule, std_log, mistakes_log


def plotly_pie_chart(input_data, legend_title):
    colors = ['gold', 'mediumturquoise', 'darkorange', 'lightgreen']
    colors = colors[:len(input_data)]
    pie_chart = go.Figure(data=[go.Pie(labels=list(input_data.keys()),
                                       values=list(input_data.values()))])
    pie_chart.update_layout(legend_title_text=legend_title)
    pie_chart.update_traces(hoverinfo='label+percent', textinfo='value', textfont_size=20,
                            marker=dict(colors=colors, line=dict(color='#000000', width=2)))
    return pie_chart


def plotly_bar_chart(input_data, legend_title=None):
    colors = ['gold', 'mediumturquoise', 'darkorange', 'lightgreen']
    colors = colors[:len(input_data)]

    bar_chart = go.Figure([go.Bar(x=list(input_data.keys()), y=list(input_data.values()))])
    # bar_chart.update_layout(legend_title_text=legend_title)
    bar_chart.update_traces(textfont_size=20,
                            marker=dict(color=colors, line=dict(color='#000000', width=2)))
    return bar_chart


def render_svg(svg):
    """Renders the given svg string. Taken from https://discuss.streamlit.io/t/display-svg/172/5
    Not safe decision, but well, c'est la vie"""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img width="400" src="data:image/svg+xml;base64,%s"/>' % b64
    return html


def get_binary_file_downloader_html(bin_file, file_label='File'):
    """Taken from https://discuss.streamlit.io/t/how-to-download-file-in-streamlit/1806/27
    Again not safe, but it's prototype"""
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">Download {file_label}</a> '
    return href


st.title('Molecules standardization with CGRtools')
text = 'This app uses [CGRtools package](https://github.com/stsouko/CGRtools) to standardize molecular data. ' \
       'It takes *SDF* or *SMILES* input files, checks their structure on valence and aromatic ring errors and ' \
       'applies standardization rules in the folowing order:  \n1. Convertation of molecule to the kekule form  \n' \
       '2. Standardization of functional groups  \n3. Hidding of explicit hydrogens  \n' \
       '4. Convertation of molecule to the aromatic (thiele) form  \n  \n' \
       'For additional information, please, refer to the [package documentation]' \
       '(https://cgrtools.readthedocs.io/tutorial/3_standardization.html)'

st.markdown(text)

input_file = st.file_uploader("Molecules in SDF or SMILES format", type=['sdf', 'smiles', 'smi'])

smiles_input = st.text_area("Input SMILES, multiple should be divided by semicolon ;", help="For example: c1ccccc1;CCO")

if st.button('Standardize molecules'):
    with st.spinner(text='Running standardization'):
        std_mistakes = {'indeces': [], 'mistakes': []}
        applied_rules_count = {'functional group': 0, 'implicify hydrogens': 0, 'aromatization': 0}
        example_flag = True
        num_std_mols = 0

        if input_file and smiles_input:
            st.error('Both file and input are given. Please, choose something one.')
            st.stop()
        elif not input_file and not smiles_input:
            st.error('Nothing was given')
            st.stop()

        with SDFWrite('standardized.sdf') as out:
            if input_file:
                file_name = input_file.name
                readed_file = input_file.getvalue().decode("utf-8")
                input_file = StringIO(readed_file)

                if file_name[-3:] == 'sdf':
                    file_length = readed_file.count('$$$$')
                    with SDFRead(input_file) as inp:
                        for n, molecule in stqdm(enumerate(inp, 1), total=file_length):
                            molecule.meta['standardized'] = False
                            new_molecule, applied_rules_count, std_mistakes = standardization_step(molecule, n,
                                                                                                   applied_rules_count,
                                                                                                   std_mistakes)
                            out.write(new_molecule)
                            if new_molecule.meta['standardized']:
                                num_std_mols += 1

                                if example_flag:
                                    example = (molecule, new_molecule,)
                                    example_flag = False

                else:
                    file_length = readed_file.count('\n') + 1
                    with SMILESRead(input_file, ignore=True) as inp:
                        for n, molecule in stqdm(enumerate(inp, 1), total=file_length):
                            molecule.meta['standardized'] = False
                            new_molecule, applied_rules_count, std_mistakes = standardization_step(molecule, n,
                                                                                                   applied_rules_count,
                                                                                                   std_mistakes)
                            out.write(new_molecule)
                            if new_molecule.meta['standardized']:
                                num_std_mols += 1

                                if example_flag:
                                    example = (molecule, new_molecule,)
                                    example_flag = False

            else:
                smiles_parser = SMILESRead.create_parser(ignore=True)
                smiles = smiles_input.strip().split(';')
                for n, smi in stqdm(enumerate(smiles, 1)):
                    molecule = smiles_parser(smi)
                    molecule.meta['standardized'] = False
                    new_molecule, applied_rules_count, std_mistakes = standardization_step(molecule, n,
                                                                                           applied_rules_count,
                                                                                           std_mistakes)
                    out.write(new_molecule)
                    if new_molecule.meta['standardized']:
                        num_std_mols += 1

                        if example_flag:
                            example = (molecule, new_molecule,)
                            example_flag = False

    if example_flag or num_std_mols == 0:
        st.warning('No molecules were standardized')
    else:
        st.subheader('Result')
        st.markdown(get_binary_file_downloader_html('standardized.sdf', 'standardized SDF file'),
                    unsafe_allow_html=True)
        st.subheader('Standardization statisitcs')

        col1, col2 = st.columns(2)
        with col1:
            if std_mistakes['mistakes']:
                df_mistakes = pd.DataFrame(std_mistakes)

                counted_mistakes = df_mistakes['mistakes'].value_counts()
                mistakes_count = {}
                all_mistakes = 0
                for i in ['invalid_aromatic_ring', 'valence_error']:
                    try:
                        mistakes_count[i] = counted_mistakes[i]
                        all_mistakes += counted_mistakes[i]
                    except KeyError:
                        mistakes_count[i] = 0
                df_mistakes.to_csv('standardization_mistakes.csv', index=False)
                mistakes_chart = plotly_pie_chart(mistakes_count, 'Errors in molecules')
                st.markdown('**Mistakes analysis**')
                st.markdown(
                    f'Your input consists of **{round(all_mistakes / file_length * 100, 2)}%** invalid molecules')
                st.plotly_chart(mistakes_chart, use_container_width=True)
                st.markdown(get_binary_file_downloader_html('standardization_mistakes.csv',
                                                            'indeces of invalid molecules'), unsafe_allow_html=True)
            else:
                st.warning('No standardization mistakes')

        with col2:
            st.markdown('**Applied rules**')
            st.markdown(f'Overall **{num_std_mols}** input molecules were standardized')

            rules_chart = plotly_bar_chart(applied_rules_count)
            st.plotly_chart(rules_chart, use_container_width=True)

        st.subheader('Example of standardization')
        col3, col4 = st.columns(2)
        with col3:
            st.write('Input molecule')
            st.write(render_svg(example[0].depict()), unsafe_allow_html=True)
        with col4:
            st.write('Standardized molecule')
            st.write(render_svg(example[1].depict()), unsafe_allow_html=True)
