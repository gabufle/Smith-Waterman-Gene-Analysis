import streamlit as st
import subprocess
import os

# --- UI SETUP ---
st.set_page_config(page_title="SW Aligner Pro", layout="centered")
st.title("🧬 Affine Smith-Waterman Aligner")
st.markdown("Powered by a custom high-performance C++ engine.")

# --- FILE UPLOADS ---
col1, col2 = st.columns(2)
with col1:
    wt_file = st.file_uploader("Upload Reference (WT) FASTA", type=["fasta", "fa", "txt"])
with col2:
    var_file = st.file_uploader("Upload Variant FASTA", type=["fasta", "fa", "txt"])

# --- THE MAGIC BUTTON ---
if st.button("Run Alignment Engine", type="primary"):
    if wt_file is None or var_file is None:
        st.warning("Please upload both sequences first!")
    else:
        with st.spinner("C++ Engine calculating affine gaps..."):
            
            # 1. Save the uploaded web files to the local hard drive temporarily
            with open("temp_wt.fasta", "wb") as f:
                f.write(wt_file.getbuffer())
            with open("temp_var.fasta", "wb") as f:
                f.write(var_file.getbuffer())
                
            # 2. Find your C++ executable (Change this if your .exe is named differently!)
            # Note: Visual Studio usually puts it in x64/Debug/
            exe_path = r"x64\Debug\Smith Waterman Alignment.exe" 
            
            if not os.path.exists(exe_path):
                st.error(f"Could not find the C++ engine at: {exe_path}")
            else:
                # 3. Command Python to run the C++ program and catch the output!
                result = subprocess.run(
                    [exe_path, "temp_wt.fasta", "temp_var.fasta"], 
                    capture_output=True, 
                    text=True
                )
                
                # 4. Paint the C++ output onto the website
                if result.returncode == 0:
                    st.success("Alignment Complete!")
                    st.code(result.stdout, language="text") # Displays that beautiful consensus line
                else:
                    st.error("Engine Error:")
                    st.code(result.stderr, language="text")
                    
            # Clean up the temporary files
            if os.path.exists("temp_wt.fasta"): os.remove("temp_wt.fasta")
            if os.path.exists("temp_var.fasta"): os.remove("temp_var.fasta")