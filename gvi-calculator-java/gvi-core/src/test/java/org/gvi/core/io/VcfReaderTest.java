package org.gvi.core.io;

import org.gvi.core.model.VariantRecord;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

class VcfReaderTest {

    @Test
    void parsesSnpAndIndelWithInfoDp() throws IOException {
        String vcf = String.join("\n",
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t100\t.\tA\tG\t60\tPASS\tDP=25",
                "chr1\t200\t.\tAT\tA\t60\tPASS\tDP=5",
                ""
        );
        VcfReader.Result result = VcfReader.read(new BufferedReader(new StringReader(vcf)), "test");
        assertThat(result.variants()).hasSize(2);
        VariantRecord snp = result.variants().get(0);
        assertThat(snp.isSnp()).isTrue();
        assertThat(snp.depth()).isEqualTo(25);
        VariantRecord indel = result.variants().get(1);
        assertThat(indel.isIndel()).isTrue();
        assertThat(indel.depth()).isEqualTo(5);
    }

    @Test
    void parsesSampleFormatDpWhenInfoDpAbsent() throws IOException {
        String vcf = String.join("\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
                "chr1\t150\t.\tC\tT\t60\tPASS\t.\tGT:DP\t0/1:42",
                ""
        );
        VcfReader.Result result = VcfReader.read(new BufferedReader(new StringReader(vcf)), "test");
        assertThat(result.variants()).hasSize(1);
        assertThat(result.variants().get(0).depth()).isEqualTo(42);
    }

    @Test
    void splitsMultiAllelicSitesIntoSeparateRecords() throws IOException {
        String vcf = String.join("\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t300\t.\tA\tG,T\t60\tPASS\tDP=10",
                ""
        );
        VcfReader.Result result = VcfReader.read(new BufferedReader(new StringReader(vcf)), "test");
        assertThat(result.variants()).hasSize(2);
    }

    @Test
    void skipsMalformedLineInsteadOfCrashing() throws IOException {
        String vcf = String.join("\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\tNOTANUMBER\t.\tA\tG\t60\tPASS\tDP=10",
                "chr1\t400\t.\tA\tG\t60\tPASS\tDP=10",
                ""
        );
        VcfReader.Result result = VcfReader.read(new BufferedReader(new StringReader(vcf)), "test");
        assertThat(result.variants()).hasSize(1);
        assertThat(result.report().getSkippedCount()).isEqualTo(1);
    }

    @Test
    void ignoresSpanningDeletionPlaceholderAllele() throws IOException {
        String vcf = String.join("\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t500\t.\tA\t*\t60\tPASS\tDP=10",
                ""
        );
        VcfReader.Result result = VcfReader.read(new BufferedReader(new StringReader(vcf)), "test");
        assertThat(result.variants()).isEmpty();
    }
}
