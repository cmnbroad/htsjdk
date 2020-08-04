package htsjdk.beta.codecs.hapref.fasta;

import htsjdk.beta.codecs.hapref.HapRefCodec;
import htsjdk.beta.codecs.hapref.HapRefDecoder;
import htsjdk.beta.codecs.hapref.HapRefEncoder;
import htsjdk.beta.plugin.bundle.Bundle;
import htsjdk.beta.plugin.registry.SignatureProbingInputStream;
import htsjdk.io.IOPath;
import htsjdk.beta.plugin.HtsDecoderOptions;
import htsjdk.beta.plugin.HtsEncoderOptions;
import htsjdk.beta.plugin.HtsCodecVersion;
import htsjdk.beta.plugin.hapref.HaploidReferenceFormat;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.utils.ValidationUtils;

import java.io.IOException;

/**
 * A FASTA codec.
 */
public class FASTACodecV1_0 extends HapRefCodec {

    public static final HtsCodecVersion VERSION_1 = new HtsCodecVersion(1, 0, 0);

    @Override
    public HtsCodecVersion getVersion() {
        return VERSION_1;
    }

    @Override
    public HaploidReferenceFormat getFileFormat() {
        return HaploidReferenceFormat.FASTA;
    }

    @Override
    public int getSignatureSize() {
        return 1;
    }

    @Override
    public boolean canDecodeSignature(final SignatureProbingInputStream rawInputStream, final String sourceName) {
        int c = rawInputStream.read();
        if (c == -1) {
            throw new RuntimeIOException(
                    String.format("Codec %s failed probing signature for resource %s", this.getDisplayName(), sourceName));
        }
        return ((char) c) == '>';  // for FASTA, this is all we have to go on...
    }

    @Override
    public boolean canDecodeURI(final IOPath ioPath) {
        return FileExtensions.FASTA.stream().anyMatch(ext-> ioPath.hasExtension(ext));
    }

   @Override
    public HapRefDecoder getDecoder(final Bundle inputBundle, final HtsDecoderOptions options) {
        ValidationUtils.validateArg(options == null, "reference reader options must be null");
        return new FASTADecoderV1_0(inputBundle);
    }

    @Override
    public HapRefEncoder getEncoder(final Bundle outputBundle, final HtsEncoderOptions options) {
        throw new IllegalStateException("Not implemented");
    }

    @Override
    public boolean runVersionUpgrade(final HtsCodecVersion sourceCodecVersion, final HtsCodecVersion targetCodecVersion) {
        throw new IllegalStateException("Not implemented");
    }
}
