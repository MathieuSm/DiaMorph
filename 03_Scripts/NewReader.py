function [image,metadata]=readISQ(filename,slicerange,window,progress)

# [image]=readISQ(filename) reads the image in the file
# [image,metadata]=readISQ(filename) also provides some metadata
# [image,metadata]=readISQ(filename,slice_nr) reads one slice at a time. slice_nr=0
# means read only metadata.
# [image]=readISQ(filename,[first_slice_nr,last_slice_nr]) reads a range of
# slices. [] means read all slices.
# [image]=readISQ(filename,slicerange,window) reads the
#   window from each slice defined by the vector
#   window=[minx,maxx,miny,maxy]. [] means read the full slice.
# [image]=readISQ(filename,slicerange,window,progress)
#   progress=1 displays progress bar (0 no bar)
#
# Johan Karlsson, AstraZeneca, 2017
# 
# See also http://www.scanco.ch/en/support/customer-login/faq-customers/faq-customers-general.html

#%%
filename = '../02_Data/C0019219_1.ISQ'
fid=open(filename, 'rb')
#h=fread(fid,128,'int') # header
h = fid.read(32) # header

#%%
metadata = {}
metadata['XDim v'] = h[12]
metadata['YDim v'] = h[13]
metadata['Zdim v'] = h[14]
metadata['XDim um'] = h[15]
metadata['YDim um'] = h[16]
metadata['ZDim um'] = h[17]
metadata['Slice Thickness um'] = h[18]
metadata['Slice Increment um'] = h[19]
metadata['Slice 1 Position um'] = h[20]
metadata['Min Data Value'] = h[21]
metadata['Max Data Value'] = h[22]
metadata['Mu Scaling'] = h[23]

metadata['Name'] = str(fid.read(40)) # header

#%%
first_slice=0
last_slice=metadata['Zdim v']
window=[0,metadata['XDim v'],0,metadata['YDim v']]
progress=0

#image=zeros(window(2)-window(1)+1,window(4)-window(3)+1,last_slice-first_slice+1)

for slice in range(first_slice, last_slice)
    fid.seek(fid,512*(h(end)+1)+metadata['XDim v']*metadata['YDim v']*2*(slice-1),'bof')
    I=permute(reshape(fread(fid,metadata.dimx_p*metadata.dimy_p,'short'),[metadata.dimx_p,metadata.dimy_p]),[2,1])
    image(:,:,slice-first_slice+1)=I(window(1):window(2),window(3):window(4))
end
if progress
    close(wb)
end
fclose(fid)
   