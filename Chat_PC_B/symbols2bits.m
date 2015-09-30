function y = symbols2bits(symbols,preamble)
% Does frame synchronization and converts symbols to bits, it is assumed
% that the symbols are in the rows of symbols (and the same for the preamble) 
% symbols - The sequence of symbols
% preamble - The preamble to be detected

[m1,n1]=size(symbols);
[m2,n2]=size(preamble);
bits=[];

for i=1:m1-m2;
    tmp=symbols(i:m2-1+i,:);
    if isequal(tmp,preamble)
        start_of_data=i+m2;
    end
end

if ~exist('start_of_data','var') %If preamble was not found, ie. start_of_data was never defined
    display('Error. Preamble not found.')
    y=0;
    return
end

for i=start_of_data:m1
    bits = [bits symbols(i,:)];
end

y = bits;

end
