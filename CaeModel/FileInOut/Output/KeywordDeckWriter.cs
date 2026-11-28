// PrePoMax - Copyright (C) 2016-2026 Matej Borovinšek
//
// Licensed under the terms defined in the LICENSE file located in the root directory of this source code.
//
// Source code: https://github.com/PrePoMax/PrePoMax
//
// Author: Matej Borovinšek
// Contributors:

using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using FileInOut.Output.Calculix;

namespace FileInOut.Output
{
    internal static class KeywordDeckWriter
    {
        // Static methods                                                                                                           
        public static void Write(string fileName, IEnumerable<CalculixKeyword> keywords, Func<string, string> format = null)
        {
            // Finish formatting before opening the destination. Opening it in place preserves its links, attributes
            // and sharing rules. The staging file is removed when its handle closes, including on a formatting error.
            string temporary = Path.Combine(Path.GetTempPath(), Path.GetRandomFileName());
            using (FileStream staged = new FileStream(temporary, FileMode.CreateNew, FileAccess.ReadWrite, FileShare.None,
                                                       4096, FileOptions.DeleteOnClose | FileOptions.SequentialScan))
            {
                using (StreamWriter writer = new StreamWriter(staged, new UTF8Encoding(false, true), 4096, true))
                {
                    WriteKeywords(writer, keywords, format);
                }
                staged.Position = 0;
                using (FileStream destination = new FileStream(fileName, FileMode.Create, FileAccess.Write, FileShare.Read))
                {
                    staged.CopyTo(destination);
                }
            }
        }
        public static void WriteKeywords(TextWriter writer, IEnumerable<CalculixKeyword> keywords, Func<string, string> format)
        {
            foreach (var keyword in keywords)
            {
                string header = keyword.GetKeywordString();
                writer.Write(format == null ? header : format(header));
                string data = keyword.GetDataString();
                writer.Write(format == null ? data : format(data));
                WriteKeywords(writer, keyword.Keywords, format);
            }
        }
    }
}
