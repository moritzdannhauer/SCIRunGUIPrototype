/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2020 Scientific Computing and Imaging Institute,
   University of Utah.

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


#ifndef SPIRE_RENDER_TEXTUREMAN_HPP
#define SPIRE_RENDER_TEXTUREMAN_HPP

#include <es-log/trace-log.h>
#include <map>
#include <set>
#include <es-cereal/CerealCore.hpp>
#include <entity-system/BaseSystem.hpp>
#include <es-systems/SystemCore.hpp>
#include <gl-platform/GLPlatform.hpp>
#include <es-acorn/Acorn.hpp>
#include <spire/scishare.h>

namespace ren {

  struct Texture;
  /// Basic texture manager. Very similar to the shader manager.
  class SCISHARE TextureMan
  {
  public:
    /// \param  numRetries  The number of retries we have to load the asset.
    ///                     Zombie promises will remain present in the system
    ///                     and a load will be re-attempted again when
    ///                     serialized and deserialized.
    explicit TextureMan(int numRetries = 2);
    ~TextureMan();

    /// Loads texture onto the given entityID.
    /// \param  core        Core base.
    /// \param  entityID    Entity ID which will receive the ren::Texture component.
    /// \param  assetName   Name of the texture to load.
    void loadTexture(spire::CerealCore& core, uint64_t entityID,
      const std::string& assetName, int32_t textureUnit,
      const std::string& uniformName);

    //create an empty texture
    ren::Texture createTexture(
      const std::string& assetName,
      GLsizei textureWidth, GLsizei textureHeight,
      GLint internalFormat, GLenum format,
      GLenum type, GLint filter);

    //create a texture from font face
    ren::Texture createTexture(
      const std::string& assetName,
      GLsizei textureWidth, GLsizei textureHeight,
      const std::vector<uint8_t>& bitmap);

    ren::Texture createTexture(
      const std::string& assetName,
      GLint internalformat,
      GLsizei width,
      GLsizei height,
      GLenum format,
      GLenum type,
      const std::vector<uint8_t>& data);

    //resize texture
    bool resizeTexture(
      ren::Texture &tex, GLsizei textureWidth,
      GLsizei textureHeight);

    /// Runs a single garbage collection cycle on the current state of the core.
    void runGCCycle(spire::ESCoreBase& core);

    /// Returns the GLID for the given assetName (texture name), if one has been
    /// generated. Returns 0 if the texture is not found in the system.
    GLuint getIDForAsset(const char* assetName) const;

    /// Retrieve the asset name from the given GLuint.
    std::string getAssetFromID(GLuint id) const;

    /// Registers TextureMan's systems. Both the garbage collector and the promise
    /// fullfillment system are registered.
    static void registerSystems(spire::Acorn& core);

    /// Obtains the garbage collectors name so that you can setup intermitent
    /// garbage collection cycles.
    static const char* getGCName();

    /// Obtains the texture promise fullfilment system's name. This system should
    /// be installed in every core. You shouldn't need to run it every frame.
    /// Maybe about every 200 MS.
    static const char* getPromiseSystemName();

  private:
    friend class TextureGarbageCollector;
    friend class TexturePromiseFulfillment;

    /// Returns false if we failed to generate the component because the asset
    /// has not been loaded yet.
    bool buildComponent(spire::CerealCore& core, uint64_t entityID,
      const std::string& assetName, int32_t textureUnit,
      const std::string& uniformName);

    /// Run garbage collection against validKeys.
    void runGCAgainstVaidIDs(const std::set<GLuint>& validKeys);

    /// Issues a request for a texture.
    void requestTexture(spire::ESCoreBase& core, const std::string& assetName,
      int32_t numRetries);

    /// Callback issued when a texture has been read from disk and is ready
    /// to load.
    void loadTextureCB(const std::string& textureName, bool error,
      size_t bytesRead, uint8_t* buffer,
      int32_t numRetries, spire::ESCoreBase& core);


    /// Called from loadTextureCB.
    void loadRawPNG(const std::string& assetName, uint8_t* buffer,
      size_t bytesRead, int numRetries,
      spire::ESCoreBase& core);

    void loadRawITX(const std::string& assetName, uint8_t* buffer,
      size_t bytesRead, int numRetries,
      spire::ESCoreBase& core);

    /// Map from GL id to asset name.
    std::map<GLuint, std::string> mGLToName;

    /// Map from asset name to GL id.
    std::map<std::string, GLuint> mNameToGL;

    /// Indicates whether new assets have been loaded but a promise fulfillment
    /// run has not been made. A GC cycle when this is true will lead to
    /// incorrect results. Hence, if a GC cycle is requested, a warning will
    /// be issued and the GC cycle will be aborted.
    bool mNewUnfulfilledAssets;

    /// Number of retries to attempt before resulting in promise failure.
    int32_t mNumRetries;
  };

} // namespace ren

#endif
