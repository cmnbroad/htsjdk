package htsjdk.variant.bcf2;

import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

/**
 */
public class BCFVCFFormatHeaderLine extends VCFFormatHeaderLine {
    final BCFVersion bcfVersion;

    public BCFVCFFormatHeaderLine(final BCFVersion bcfVersion, String line, VCFHeaderVersion version)
    {
        super(line, version);
        if (bcfVersion.getMinorVersion() < BCFVersion.BCF_MINOR_VERSION_2) {
            throw new IllegalStateException(
                    String.format("BCFVCFFormatHeaderLine only valid for BCF versions > " + bcfVersion));
        }
        this.bcfVersion = bcfVersion;
    }

    @Override
    protected boolean getIDXFieldAllowed() {
        return true;
    }

}
