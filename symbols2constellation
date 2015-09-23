function y = symbols2constellation(symbols, const)
% Takes a vector of symbols and converts it into the constellation const

[m,n]=size(symbols);

if strcmp(const,'4QAM')
    if n==2
        % Create a 4-QAM constellation in which symbol 00 corresponds to positive I and positive Q (++), 10 to -+, 11 to --, 01 to +-
        s=[(1 + 1i) (1 - 1i) (-1 +1i) (-1 - 1i)]/sqrt(2);
        symbol_index = bi2de(symbols, 'left-msb')'+1; % Symbols to symbol index
        y = s(symbol_index); % Put the symbols in the constellation
    else
        display('Error. Symbols have the wrong format for 4QAM.');
        y=0;
        return;
    end
else
    display('Error. Invalid constellation.');
    y=0;
    return;
end

end
