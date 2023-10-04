from Bio import Align
from Bio.Seq import Seq
import requests
import pandas as pd
import numpy as np
import re
import streamlit as st 
# Creating sample sequences
with st.form(key='parameters'):
    reference=st.text_input('Enter reference transcript:','NM_001429.3',help='Example: NM_001429.3')
    variant = st.text_input('Enter variant', 'NM_001429.3(EP300_v001):c.4810_4886dup',help='Example: NM_001429.3(EP300_v001):c.4810_4886dup')  
    submit_button = st.form_submit_button(label='Submit')
if submit_button:

    mutalyzer_variant_seq=requests.get("https://mutalyzer.nl/api/mutate/"+variant)
    seq_var=mutalyzer_variant_seq.json()['sequence']['seq']

    mutalyzer_reference_seq = requests.get('https://mutalyzer.nl/api/reference_model/?reference_id='+reference+'&siblings=false&ancestors=false&descendants=false')
    seq_ref=mutalyzer_reference_seq.json()['sequence']['seq']

    # Finding similarities
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -2
    alignments = aligner.align(seq_ref, seq_var)
    #seqshow1= [a for a in str(alignments[0].sequences[0])] 
    #seqshow2= [a for a in str(alignments[0].sequences[1])]
    ref_a=re.search('-',alignments[0]._get_row(0))
    var_a=re.search('-',alignments[0]._get_row(1))
    if ref_a != None and var_a != None:
        rr=ref_a.span()+var_a.span()
        min_val, max_val = np.min(rr), np.max(rr)
        max_val=max_val+np.max((len(re.findall('-',alignments[0]._get_row(0))),len(re.findall('-',alignments[0]._get_row(1)))))
    elif ref_a != None:
        rr=ref_a.span()
        min_val, max_val = np.min(rr), np.max(rr)
        max_val=max_val+len(re.findall('-',alignments[0]._get_row(0)))
    else:
        rr=var_a.span()   
        min_val, max_val = np.min(rr), np.max(rr)
        max_val=max_val+len(re.findall('-',alignments[0]._get_row(1))) 
    array1=[]
    rc_array1=[]
    c_array1=[]
    r_array1=[]
    array2=[]
    rc_array2=[]
    c_array2=[]
    r_array2=[]    
    
    def rev(it):
        "Reverses an interable and returns it as a string"
        return ''.join(reversed(it))
    
    for ch2 in alignments[0]._get_row(1)[min_val-5:max_val+5]:
        array2.append(ch2)
    for ch1 in alignments[0]._get_row(0)[min_val-5:max_val+5]:
        array1.append(ch1)    
    for rc_ch2 in str(Seq(alignments[0]._get_row(1)[min_val-5:max_val+5]).reverse_complement()):
        rc_array2.append(rc_ch2)
    for rc_ch1 in str(Seq(alignments[0]._get_row(0)[min_val-5:max_val+5]).reverse_complement()):
        rc_array1.append(rc_ch1)
    for c_ch2 in str(Seq(alignments[0]._get_row(1)[min_val-5:max_val+5]).complement()):
        c_array2.append(c_ch2)
    for c_ch1 in str(Seq(alignments[0]._get_row(0)[min_val-5:max_val+5]).complement()):
        c_array1.append(c_ch1)
    for r_ch2 in rev(str(Seq(alignments[0]._get_row(1)[min_val-5:max_val+5]))):
        r_array2.append(r_ch2)
    for r_ch1 in rev(str(Seq(alignments[0]._get_row(0)[min_val-5:max_val+5]))):
        r_array1.append(r_ch1)

    towrite1=pd.concat((pd.DataFrame(array1).T,pd.DataFrame(array2).T))
    towrite2=pd.concat((pd.DataFrame(rc_array1).T,pd.DataFrame(rc_array2).T))
    towrite3=pd.concat((pd.DataFrame(c_array1).T,pd.DataFrame(c_array2).T))
    st.write('mutalyzer sequence')
    st.write(towrite1)
    chart_data1 = towrite1.apply(lambda x: x.iloc[0]==x.iloc[1])
    st.bar_chart(chart_data1.T,width=towrite1.shape[1])

    st.write('reverse complement')
    st.write(towrite2)
    chart_data2 = towrite2.apply(lambda x: x.iloc[0]==x.iloc[1])
    st.bar_chart(chart_data2.T,width=towrite2.shape[1])
    st.write('complement')
    st.write(towrite3)
    chart_data3 = towrite3.apply(lambda x: x.iloc[0]==x.iloc[1])
    st.bar_chart(chart_data3.T,width=towrite3.shape[1])

    st.write('reverse')
    towrite4=pd.concat((pd.DataFrame(r_array1).T,pd.DataFrame(r_array2).T))
    st.write(towrite4)
    chart_data4 = towrite4.apply(lambda x: x.iloc[0]==x.iloc[1])
    st.bar_chart(chart_data4.T,width=towrite4.shape[1])
