package htsjdk.variant.bcf2;

import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.vcf.*;

/**
 */
public class BCFVCFCodec extends VCFCodec {

    final BCFVersion bcfVersion;

    BCFVCFCodec(final BCFVersion version) {
        if (version == null) {
            throw new TribbleException("Null version in BCFVCFCodec");
        }
        bcfVersion = version;
    }

    /**
     * get the VCFHeaderVersion from this BCF's VCF header
     * @return VCFHeaderVersion
     */
    public VCFHeaderVersion getVCFVersion() { return version; }

    @Override
    public Object readActualHeader(final LineIterator lineIterator) {
        return super.readActualHeader(lineIterator);
    }

    /**
     * Allow subclasses to override this in order to provide custom FILTER line handling.
     *
     * @param line
     * @param version
     * @return
     */
    @Override
    public VCFFilterHeaderLine getFilterHeaderLineHandler(final String line, final VCFHeaderVersion version) {
        return bcfVersion.getMinorVersion() >= BCFVersion.BCF_MINOR_VERSION_2 ?
            new BCFVCFFilterHeaderLine(bcfVersion, line, version) :
            new VCFFilterHeaderLine(line, version);
    }

    /**
     * Allow subclasses to override this in order to provide custom INFO line handling.
     *
     * @param line
     * @param version
     * @return
     */
    @Override
    public VCFInfoHeaderLine getInfoHeaderLineHandler(final String line, final VCFHeaderVersion version) {
        return bcfVersion.getMinorVersion() >= BCFVersion.BCF_MINOR_VERSION_2 ?
            new BCFVCFInfoHeaderLine(bcfVersion, line, version) :
            new VCFInfoHeaderLine(line, version);
    }

    /**
     * Allow subclasses to override this in order to provide custom INFO line handling.
     *
     * @param line
     * @param version
     * @return
     */
    @Override
    public VCFFormatHeaderLine getFormatHeaderLineHandler(final String line, final VCFHeaderVersion version) {
        return bcfVersion.getMinorVersion() >= BCFVersion.BCF_MINOR_VERSION_2 ?
            new BCFVCFFormatHeaderLine(bcfVersion, line, version) :
            new VCFFormatHeaderLine(line, version);
    }

}
