function suffix = featureLabelChecksum(text)
%FEATURELABELCHECKSUM CRC-32/ISO-HDLC of UTF-8, only for bounded display labels.
% This eight-hex checksum is not collision-free, cryptographic, or a seed mixer.
% It does not identify an experiment: callers retain the full label and use the
% existing unique directory allocator. Integer operations need no JVM or tools.
    crc = uint32(hex2dec('FFFFFFFF'));
    polynomial = uint32(hex2dec('EDB88320'));
    bytes = unicode2native(text, 'UTF-8');
    for i_byte = 1:numel(bytes)
        crc = bitxor(crc, uint32(bytes(i_byte)));
        for i_bit = 1:8
            if bitand(crc, uint32(1))
                crc = bitxor(bitshift(crc, -1), polynomial);
            else
                crc = bitshift(crc, -1);
            end
        end
    end
    suffix = lower(dec2hex(bitcmp(crc), 8));
end
