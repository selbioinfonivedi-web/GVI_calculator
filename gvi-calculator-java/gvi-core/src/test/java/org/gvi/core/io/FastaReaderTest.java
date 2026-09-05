package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.NucleotideSequence;
import org.junit.jupiter.api.Test;

import java.io.StringReader;
import java.time.LocalDate;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

class FastaReaderTest {

    @Test
    void parsesSimpleMultiFasta() {
        String fasta = ">seq1 description here\nACGT\nACGT\n>seq2\nTTTT\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        assertThat(result.sequences()).hasSize(2);
        assertThat(result.sequences().get(0).getId()).isEqualTo("seq1");
        assertThat(result.sequences().get(0).getSequence()).isEqualTo("ACGTACGT");
        assertThat(result.sequences().get(1).getSequence()).isEqualTo("TTTT");
        assertThat(result.report().getSkippedCount()).isZero();
    }

    @Test
    void skipsDuplicateIdsWithWarningInsteadOfCrashing() {
        String fasta = ">seq1\nACGT\n>seq1\nGGGG\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        assertThat(result.sequences()).hasSize(1);
        assertThat(result.sequences().get(0).getSequence()).isEqualTo("ACGT");
        assertThat(result.report().getSkippedCount()).isEqualTo(1);
        assertThat(result.report().getWarnings().get(0).reason()).contains("Duplicate");
    }

    @Test
    void skipsEmptySequenceRecord() {
        String fasta = ">empty\n>seq2\nACGT\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        assertThat(result.sequences()).hasSize(1);
        assertThat(result.report().getSkippedCount()).isEqualTo(1);
    }

    @Test
    void skipsSequenceWithInvalidCharacters() {
        String fasta = ">bad\nACGT123XYZ!!\n>good\nACGT\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        assertThat(result.sequences()).hasSize(1);
        assertThat(result.sequences().get(0).getId()).isEqualTo("good");
        assertThat(result.report().getWarnings().get(0).reason()).contains("IUPAC");
    }

    @Test
    void throwsTypedExceptionOnEmptyFile() {
        assertThatThrownBy(() -> FastaReader.read(new StringReader(""), "empty-file"))
                .isInstanceOf(GviInputException.class);
    }

    @Test
    void handlesAmbiguityCodesAndGaps() {
        String fasta = ">withAmbig\nACGTNRYWSKMBDHV--\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        assertThat(result.sequences()).hasSize(1);
        List<NucleotideSequence> seqs = result.sequences();
        assertThat(seqs.get(0).length()).isEqualTo(17);
    }

    @Test
    void extractsCollectionDateFromGisaidStyleHeaderAutomatically() {
        String fasta = ">hCoV-19/USA/CA-1/2020|2020-03-15|EPI_ISL_12345\nACGT\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        assertThat(result.sequences().get(0).getCollectionDate()).contains(LocalDate.of(2020, 3, 15));
    }

    @Test
    void sequenceWithNoDateInHeaderHasEmptyCollectionDate() {
        String fasta = ">plain_sample_id\nACGT\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        assertThat(result.sequences().get(0).getCollectionDate()).isEmpty();
    }

    @Test
    void extractsLocationAndHostAlongsideDateFromTheAccessionLocationHostDateConvention() {
        String fasta = ">KP821387.1|France|domestic_sheep|2001\nACGT\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        NucleotideSequence seq = result.sequences().get(0);
        // Midpoint of the year, not January 1st -- see TemporalUtil#midpointOfYear.
        assertThat(seq.getCollectionDate()).contains(org.gvi.core.util.TemporalUtil.midpointOfYear(2001));
        assertThat(seq.getLocation()).contains("France");
        assertThat(seq.getHost()).contains("domestic_sheep");
    }

    @Test
    void doesNotExtractLocationOrHostFromAGisaidStyleHeaderThatDoesNotMatchTheFourFieldConvention() {
        String fasta = ">hCoV-19/USA/CA-1/2020|2020-03-15|EPI_ISL_12345\nACGT\n";
        FastaReader.Result result = FastaReader.read(new StringReader(fasta), "test");
        NucleotideSequence seq = result.sequences().get(0);
        assertThat(seq.getCollectionDate()).contains(LocalDate.of(2020, 3, 15)); // date still extracted
        assertThat(seq.getLocation()).isEmpty(); // but location/host are not guessed at
        assertThat(seq.getHost()).isEmpty();
    }
}
