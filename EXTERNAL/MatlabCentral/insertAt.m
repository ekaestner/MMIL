function arrOut = insertAt(arr,val,index)

arrOut = arr;
for iAS = 1:numel(index)
    assert( index(iAS)<= numel(arr)+1);
    assert( index(iAS)>=1);
    if index(iAS) == numel(arrOut)+1
        arrOut = [arrOut val(iAS)];
    else
        arrOut = [arrOut(1:index(iAS)-1+(iAS-1)) val(iAS) arrOut(index(iAS)+(iAS-1):end)];
    end
end

end
