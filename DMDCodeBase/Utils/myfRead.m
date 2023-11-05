function readNum = myfRead(port,n)
% faster read function than matlab's as it does not wait for a nextline
if (port.BytesAvailable > 0)
    readVal = fread(port, port.BytesAvailable);
    readNum = char(readVal(1:n))';
else
    readNum = -1;
end
end