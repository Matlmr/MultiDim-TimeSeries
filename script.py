from serialTimer import SerialTimer

st = SerialTimer(aoi_file='/home/mlamarre/Documents/FieldTimeSeries/outline.csv', out_file='/home/mlamarre/Documents/FieldTimeSeries/dataframes/austria.csv',dimensions=['BS.VH','BS.VV','coh.VH','coh.VV','ha_alpha.Alpha','ha_alpha.Anisotropy','ha_alpha.Entropy'], quick_check=False)

st.routine()