package htsjdk.variant.bcf2;

import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 */
public class BCFVCFInfoHeaderLine extends VCFInfoHeaderLine{

    final BCFVersion bcfVersion;

    BCFVCFInfoHeaderLine(final BCFVersion bcfVersion, final String line, final VCFHeaderVersion version) {
        super(line, version);
        if (bcfVersion.getMinorVersion() < BCFVersion.BCF_MINOR_VERSION_2) {
            throw new IllegalStateException(
                    String.format("BCFVCFInfoHeaderLine only valid for BCF versions > " + bcfVersion));
        }
        this.bcfVersion = bcfVersion;
    }

    @Override
    protected boolean getIDXFieldAllowed() {
        return true;
    }

}
