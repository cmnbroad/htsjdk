package htsjdk.variant.bcf2;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.util.Arrays;

public class BCFVCFFilterHeaderLine extends VCFFilterHeaderLine {
    final BCFVersion bcfVersion;

    BCFVCFFilterHeaderLine(final BCFVersion bcfVersion, final String line, final VCFHeaderVersion version) {
        // include IDX
        super(line, version, Arrays.asList("ID", "Description", BCF2Codec.IDXField));
        if (bcfVersion.getMinorVersion() < BCFVersion.BCF_MINOR_VERSION_2) {
            throw new IllegalStateException(
                    String.format("BCFVCFFormatHeaderLine only valid for BCF versions > " + bcfVersion));
        }
        this.bcfVersion = bcfVersion;
    }

}
