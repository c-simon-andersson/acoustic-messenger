function y = samples2symbols(iWaveform, qWaveform, const)
% Converts the iWaveform and qWaveform to symbols
% iWaveform and qWaveform contains an amount of samples
% const is the constellation used to convert to symbols

s=[];

if length(iWaveform)~=length(qWaveform)
    display('Error. iWaveform and qWaveform are not equal in length.')
    y=0;
    return;
else
    for i=1:length(iWaveform)
        s=[s; iWaveform(i)+1i*qWaveform(i)];
    end
    if strcmp(const,'4QAM')
        symbols=[];
        for j=1:length(s)
            if real(s(j))>=0 && imag(s(j))>=0
                symbols=[symbols; 0 0];
            elseif real(s(j))>=0 && imag(s(j))<0
                symbols=[symbols; 0 1];
            elseif real(s(j))<0 && imag(s(j))>=0
                symbols=[symbols; 1 0];
            elseif real(s(j))<0 && imag(s(j))<0
                symbols=[symbols; 1 1];
            else
                Display('Error. Something went wrong when mapping to symbols.');
            end
        end
    else
        display('Error. Invalid constellation.');
        y=0;
        return;
    end
end

y=symbols;

end
