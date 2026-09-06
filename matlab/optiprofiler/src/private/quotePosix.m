function quoted = quotePosix(value)
%QUOTEPOSIX Single-quote one literal path, including embedded apostrophes.
    quoted = char("'" + replace(string(value), "'", "'\''") + "'");
end
