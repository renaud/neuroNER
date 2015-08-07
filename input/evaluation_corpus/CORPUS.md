# neuroNER evaluation corpus

This corpus was created to evaluate neuroNER's performance (precision, recall). It consists of 97 pdfs (see `corpus_pmids.txt`) cited in an upcoming HBP article. 

## Preprocessing: bluima pipeline script 

Using `bluima_20150629`

    cr: ch.epfl.bbp.uima.pdf.cr.PdfCollectionReader
     inputDirectory: /Volumes/HDD2/ren_data/data_hdd/_papers_etc/bbp_people/20130701_pdf_references_from_column_paper/
     expandAbbrevs__java: true
    ae_java: ch.epfl.bbp.uima.ae.OpenNlpHelper.getSentenceSplitter();
    ae: SentenceDumpAnnotator
    ae: StatsAnnotatorPlus
     printEvery__java: 1

Results in `bluima_20150629/sentence.txt`, saved here as `corpus_raw_sentences.txt`

## Preprocessing: manual steps

* remove sentences not starting with a pmid (`^(?!\d{4})`)
* removed sentences not ending with a dot (`[^\.]$`)
* removed short sentences `^.{0,60}$`
* kept only sentences with `^.*(neuron|cell)s?.*$`
* removed lines with citation dates `\(19\d\d[abcd]?\)` and `\(20[01]\d[abcd]?\)`
* shuffled sentences

